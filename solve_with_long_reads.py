#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 12:37:18 2021

@author: rfaure
"""

import re #to find all numbers in a mixed number/letters string (such as 31M1D4M), to split on several characters (<> in longReads_interactionMatrix)
from collections import Counter #to count the elements in a list quickly
from scipy import sparse #to handle sparse matrices
from determine_multiplicity import determine_multiplicity

import segment

#Master function of the file
#Input : initial gfa (as a list of segments), a GAF file with long reads mapped to the segments, names (which is an index numbering the contig), multiplicity : the pre-computed ploidy of each contig (as numbered in names)
#Output : new gfa (as a list of segments) corrected with long reads, and modified copiesnumber (taking into account contigs that have been duplicated)
def bridge_with_long_reads(segments, names, copiesnumber, gafFile):
    
    ##There are two phases : first build consensus with approximate haploid contigs. Detect the incoherence, obtain a reliable haploid list and do that all over again
    
    #first phase : determine all contigs that have strictly one link at each end (they are approximately haploid), and sort them by length
    
    refHaploidy = determine_multiplicity(segments, names)
    haploidContigs = []
    for s in segments :
        
        #to be deemed haploid, a segment must have at mosetone connection at each of its end plus be equally or less covered thant neighboring contigs
        links = s.links
        if round(s.depths[0]/refHaploidy) == 1 :
            haploidContigs.append(s)
        elif len(links[0]) == 1 and len(links[1]) == 1 :# and s.depths[0] < 1.5*links[0][0].depths[0] and s.depths[0] < 1.5*links[1][0].depths[0] : 
            haploidContigs.append(s)
        elif len(links[0]) == 1 and len(links[1]) == 0  : #and s.depths[0] < 1.5*links[0][0].depths[0] : 
            haploidContigs.append(s)
        elif len(links[0]) == 0 and len(links[1]) == 1 : # and s.depths[0] < 1.5*links[1][0].depths[0] : 
            haploidContigs.append(s)
    
    haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in haploidContigs :
        haploidContigsNames[s.names[0]] = index
        index += 1

    #inventoriate all bridges in the list bridges : sequence of contigs found in the GAF containing at least one contig of multiplicity 1
    bridges = [[[],[]] for i in range(len(haploidContigs))] #bridges is a list inventoring at index haploidCOntigsNames[seg.names[0]] all the links left and right of the contig, supported by the gaf
    inventoriate_bridges(bridges, haploidContigsNames, gafFile, 0.9, 0.7) 
    
    #now, from all the bridges, build consensus bridges
    consensus_bridges = [['',''] for i in range(len(haploidContigs))] #consensus bridge is essentially the same as bridges, except there is only one bridge left at each side for each contig
    supported_links = sparse.dok_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file
    build_consensus_bridges(consensus_bridges, bridges, names, supported_links, haploidContigs)
    
    #consensus bridges overlap two by two (e.g. >2>3>4 right of 1 overlaps with <3<2<1 left of 4), so merge them to have a set of non-overlapping consensus bridges
    non_overlapping_bridges = [['',''] for i in range(len(haploidContigs))] 
    reliable_haploid_contigs, reliable_haploid_contigsNames = merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs)

    # print(reliable_haploid_contigsNames)
    # print("Phase 1 over...\n\n")
    #the first phase is over, now we have the reliable_haploid_contigs : let's re-do everything with these more reliable haploid contigs
    bridges = [[[],[]] for i in range(len(reliable_haploid_contigs))] #bridges is a list inventoring at index haploidCOntigsNames[seg.names[0]] all the links left and right of the contig, supported by the gaf
    inventoriate_bridges(bridges, reliable_haploid_contigsNames, gafFile, 0.9, 0.7)
        
    supported_links = sparse.dok_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file
    consensus_bridges = [['',''] for i in range(len(reliable_haploid_contigs))] #consensus bridge is essentially the same as bridges, except there is only one bridge left at each side for each contig
    build_consensus_bridges(consensus_bridges, bridges, names, supported_links, reliable_haploid_contigs)

    non_overlapping_bridges = [['',''] for i in range(len(reliable_haploid_contigs))] 
    merge_bridges(non_overlapping_bridges, consensus_bridges, reliable_haploid_contigsNames, reliable_haploid_contigs)
    #print(non_overlapping_bridges[reliable_haploid_contigsNames['134']])
    
    #second phase is over
    
    #now actually unzip the graph using the instructions in non_overlapping_bridges
    unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, reliable_haploid_contigs, reliable_haploid_contigsNames, names, supported_links)
        
    #print(non_overlapping_bridges)
    
    return segments, copiesnumber

#input : the name of the gaf file, as well as two parameters to know which alignment to keep (similarity threshold and proportion of the lenght of the query that mathched)
#output : the completed bridges list, with for each haploid contig a list of what was found left and right of the contig
def inventoriate_bridges(bridges, haploidContigsNames, gafFile, similarity_threshold, whole_mapping_threshold) :
    
    gaf = open(gafFile, 'r')
    
    for line in gaf :
        ls = line.split('\t')
        path = ls[5] # in GAF format, the 6th column is the path on which the read matched
        

        if ls[5].count('>') + ls[5].count('<') > 1 :
            
            

            if (not 'id:f' in ls[-2]) or (float(ls[-2].split(':')[-1]) > similarity_threshold) :
                
                
                if float(ls[1]) / (float(ls[3])-float(ls[2])) > whole_mapping_threshold :
                    
                    
                    contigs = re.split('[><]' , ls[5])
                    orientations = "".join(re.findall("[<>]", ls[5]))
                    del contigs[0] #because the first element is always ''

                    for c, contig in enumerate(contigs) :
                        
                        if contig in haploidContigsNames :
                            
                            
                            if orientations[c] == ">" :
                                r = 0
                                #first look at what contigs are left of the contig of interest
                                bridges[haploidContigsNames[contig]][1] +=  [""]
                                for c2 in range(c+1, len(contigs)) :
                                    
                                    bridges[haploidContigsNames[contig]][1][-1] += orientations[c2] + contigs[c2]
                                    
                                    if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                                        break
                                    
                                #then look at what's left of the contig of interest (so mirror the orientations)
                                bridges[haploidContigsNames[contig]][0] +=  [""]
                                for c2 in range(c-1 , -1, -1) :
                                    
                                    if orientations[c2] == '>' :
                                        bridges[haploidContigsNames[contig]][0][-1] += '<' + contigs[c2]
                                    else :
                                        bridges[haploidContigsNames[contig]][0][-1] += '>' + contigs[c2]
                                        
                                    if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                                        break
                                        
                            else :
                                
                                #first look at what contigs are left of the contig of interest
                                bridges[haploidContigsNames[contig]][0] +=  [""]
                                for c2 in range(c+1, len(contigs)) :
                                    bridges[haploidContigsNames[contig]][0][-1] += orientations[c2] + contigs[c2]
                                    if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                                        break
                                    
                                #then look at what's left of the contig of interest (so mirror the orientations)
                                bridges[haploidContigsNames[contig]][1] +=  [""]
                                for c2 in range(c-1 , -1, -1 ) :
                                    
                                    if orientations[c2] == '>' :
                                        bridges[haploidContigsNames[contig]][1][-1] += '<' + contigs[c2]
                                    else :
                                        bridges[haploidContigsNames[contig]][1][-1] += '>' + contigs[c2]
                                        
                                    if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                                        break
                                    
#input : list of bridges for each haploid contig
#output : completed consensus_bridges, where there is max one bridge at each end of contig
def build_consensus_bridges(consensus_bridges, bridges, names, supported_links, haploidContigs):
    
    #delete all the empty strings in bridges, they will interfere in consensus building
    bridges =  [ [ [i[j][k] for k in range(len(i[j])) if i[j][k] != "" ] for j in range(2)] for i in bridges]
    allcontigs = [ [ [ re.split('[><]' , i[j][k])[1:] for k in range(len(i[j])) ] for j in range(2)] for i in bridges]
    allorientations = [ [ [ "".join(re.findall("[<>]", i[j][k])) for k in range(len(i[j])) ] for j in range(2)] for i in bridges]
    
    #print(allcontigs[0][1])
    
    for c in range(len(bridges)) :
        
        for end in range(2) : #look left and right of each contig
        
            kept_reads = [i for i in range(len(bridges[c][end]))] #list keeping in mind all the reads still in the consensus
            pos = 0
            consensus = ''
            previousContig = names[haploidContigs[c].names[0]] # variable containing the index of the last contig of the consensus (useful to fill supported_links)
            previousEnd = end #variable containing the free end of the last contig of the consensus (useful to fill supported_links)
            
            while len(kept_reads) > 0 :
                
                candidate2 = [allcontigs[c][end][i][pos] for i in kept_reads]
                candidate1 = [allorientations[c][end][i][pos] for i in kept_reads]
                
                cons1 = Counter(candidate1).most_common(1)[0]
                cons2 = Counter(candidate2).most_common(1)[0]
                                
                if cons1[1] > 0.8*len(kept_reads) and cons2[1] > 0.8*len(kept_reads)  : #consensus is defined there as a 80% of the propositions agreeing
                
                    if cons1[0] == '*' :
                        break
                    
                    consensus_bridges[c][end] = consensus_bridges[c][end] + cons1[0] + cons2[0]
                    
                    #inventoriate this bridge in supported_links
                    
                    #print('Adding from ', consensus_bridges[c][end], " c: ", c)
                    current_end = 0
                    if cons1[0] == "<" :
                        current_end=1
                    supported_links[2*names[cons2[0]]+current_end, 2*previousContig+previousEnd] += 0.5 #since each link is present twice in links, only add 0.5
                    supported_links[2*previousContig+previousEnd, 2*names[cons2[0]]+current_end] += 0.5
                    previousContig = names[cons2[0]]
                    previousEnd = 1-current_end
                    
                    #update the list of reads still there to build the consensu
                    new_kept = []
                    for r in kept_reads:
                        if len(allcontigs[c][end][r]) > pos+1 :
                            if allorientations[c][end][r][pos] == cons1[0] and allcontigs[c][end][r][pos] == cons2[0] :
                                new_kept += [r]
                        
                    kept_reads = new_kept
               
                else : #if there is no consensus let's stop
                    break
                
                pos += 1
                                
#input: list of all consensus bridges for all haploid contigs, plus a list of the haploid contigs
#output: a list of non-overlapping bridges left and right of each contig, with only full bridges (connecting 2 haploid contigs)
def merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs) :
    
    not_actually_haploid = [] #a list of not actually haploid contigs among the haploidContigs
    for c in range(len(consensus_bridges)) :
        
        for end in range(2) :
            
            if consensus_bridges[c][end] != '' :
            
                contigs = re.split('[><]' , consensus_bridges[c][end])
                orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
                del contigs[0] #because the first element is always ''
                
                if contigs[-1] in haploidContigsNames : #in case there is a full bridge
                    #check if the two bridges are coherent
                    coherent = False
                    otherEnd = 0
                    if orientations[-1] == "<":
                        otherEnd = 1
                        
                        # print ("here is contigs [0] : ", haploidContigs[c].names[0], " -- ", re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[-1]]][otherEnd]))

                    if haploidContigs[c].names[0] == re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[-1]]][otherEnd])[-1] :
                        coherent = True
                            
                    if not coherent :
                        if c < haploidContigsNames[contigs[-1]] :
                            #print("two incoherent bridges between ", haploidContigs[c].names[0] , " and ", haploidContigs[haploidContigsNames[contigs[-1]]].names[0], " : ", consensus_bridges[c][end], " vs ", consensus_bridges[haploidContigsNames[contigs[-1]]][otherEnd])
                            not_actually_haploid += [haploidContigsNames[contigs[-1]]]
                        consensus_bridges[c][end] = ''
                        
                    else :
                        non_overlapping_bridges[c][end] = consensus_bridges[c][end]
                        consensus_bridges[c][end] = '' 
                
                else : #a half-bridge cannot be unambiguously extacted from the graph
                    non_overlapping_bridges[c][end] = ''
    
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    not_actually_haploid = list(set(not_actually_haploid))
    not_actually_haploid.sort()
        
    index = 0
    indexNot = 0
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(not_actually_haploid) and i == not_actually_haploid[indexNot] :
            indexNot += 1
        else :
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
        
    return reliable_haploid_contigs, reliable_haploid_contigsNames #return the updated list of haploid contigs
          
                
#input : a list of segments and the non_overlapping_bridges
#output: a list of segment where all the bridges have been built          
def unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, haploidContigs, haploidContigsNames, names, supported_links) :

    #compute the minimum multiplicity of each contig, so that there is enough contig to build all bridges, even when the depth suggests otherwise
    minimum_multiplicity = [0 for i in range(len(names))] #list that will contain the minimum multiplicity of each contigs, i.e. the number of different bridges in which it is found
    for s, seg in enumerate(segments) :
        minLeft = 0
        minRight = 0
        
        for n, neighbor in enumerate(seg.links[0]) :
            minLeft += supported_links[2*names[seg.names[0]], 2*names[neighbor.names[0]] + seg.otherEndOfLinks[0][n]]
            
        for n, neighbor in enumerate(seg.links[1]) :
            minRight += supported_links[2*names[seg.names[0]]+1, 2*names[neighbor.names[0]] + seg.otherEndOfLinks[1][n]]
        
        minimum_multiplicity[names[seg.names[0]]] = max(minLeft, minRight)
    
    l = len(segments)
    for se in range(l) :
        
        s = segments[se]
        
        if s.names[0] in haploidContigsNames :
            
            for end in range(2) :
                
                
                if len(non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]) > 0 : #means there is a bridge to be built there
                    
                    contigs = re.split('[><]' , non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                    orientations = "<>"[end] + "".join(re.findall("[<>]", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]))
                    del contigs[0] #because the first element is always ''
                    contigs = [s.names[0]] + contigs
                    contigsID = [segments[names[c]].ID for c in contigs]
                    newContigsIndices = [names[contigs[0]]] #list keeping track of the indices of the contigs on the path
                    oldContigsIndices = [names[i] for i in contigs]
                    
                    haploidCoverage = (segments[names[contigs[0]]].depths[0]+segments[names[contigs[-1]]].depths[0])/2 #computation of the haploid coverage at this point in the assembly

                    #if '24' in contigs :
                    
                    #print("let's duplicate ", contigs, " ", orientations)
                    nextEnd = 0
                    if orientations[1] == '<' :
                        nextEnd = 1
                    CIGAR = segments[names[contigs[0]]].CIGARs[end][segment.find_this_link(segments[oldContigsIndices[1]], nextEnd, segments[names[contigs[0]]].links[end], segments[names[contigs[0]]].otherEndOfLinks[end])]
                    nextCIGAR = '';
                    
                    for c in range(1, len(contigs)) :
                        
                        contig = segments[names[contigs[c]]]
                            
                        multiplicity = max([1, round(contig.depths[0]/haploidCoverage), minimum_multiplicity[names[contigs[c]]]]) #the multiplicity is inferred from the coverage of the contig compared to the coverage of the two haploid contigs at the extremity of the bridge
                        multiplicity = minimum_multiplicity[names[contigs[c]]]
                        minimum_multiplicity[names[contigs[c]]] -= 1
                        neighborMultiplicity = max([1, round(segments[oldContigsIndices[c-1]].depths[0]/haploidCoverage), minimum_multiplicity[names[contigs[c-1]]]])
                        neighborMultiplicity = minimum_multiplicity[names[contigs[c-1]]]
                        
                        
                        end1, end0 = 1, 1
                        if orientations[c] == '>' :
                            end1 = 0
                        if orientations[c-1] == '<' :
                            end0 = 0
                            
                        if c < len(contigs)-1 :
                            nextEnd = 0
                            if orientations[c+1] == '<' :
                                nextEnd = 1
                            nextCIGAR = contig.CIGARs[1-end1][segment.find_this_link(segments[oldContigsIndices[c+1]], nextEnd, contig.links[1-end1], contig.otherEndOfLinks[1-end1])]
                            
                        #print("Link supported with strength ", supported_links[names[contigs[c]*2+end1 , names[contigs[c-1]]*2+end0])
                                                    
                        if multiplicity > 1 and c < len(contigs)-1: #if multiplicity>1, the contig should be duplicated 
                            newSegment = segment.Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i/multiplicity for i in contig.depths])
                            segments.append(newSegment)
                            newContigsIndices += [len(segments) - 1]
                            
                            for n in newSegment.names :
                                copiesnumber[n] += 1

                            #add the link to form the new bridge
                            segment.add_link(segments[-1], end1, segments[newContigsIndices[c-1]], end0, CIGAR)
                            
                            #delete the old link if and only if it was only supported by one path
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            if supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] == 0 and segment.find_this_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]].links[end0], segments[oldContigsIndices[c-1]].otherEndOfLinks[end0], warning=False) != -1:
                                #print("remove links left of ", segments[oldContigsIndices[c]].names, " : ", segments[oldContigsIndices[c-1]].names)
                                segment.delete_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]], end0, warning = True) #though it is very hard to be sure, it is not impossible at that point that we actually delete a link that is present in another bridge, so don't warn
                            #since contig has been duplicated, lower its depth
                            contig.divide_depths(multiplicity/(multiplicity-1))
                                                        
                        else :
                            
                            # if the contig has multiple incoming reads left, get rid of all of them except (the good one is reestablished later). Do the same at the other end of the contig
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            if len(contig.links[end1]) > 0 :
                                for n, neighbor in enumerate(contig.links[end1]) :
                                    if neighbor.ID != segments[newContigsIndices[c-1]].ID:
                                        if supported_links[names[contigs[c]]*2+end1 , names[neighbor.names[0]]*2+contig.otherEndOfLinks[end1][n]] == 0 :
                                        #print("deleting links left of ", contig.names, " : ", neighbor.names, " :: ", [i.names for i in contig.links[end1]])
                                            segment.delete_link(contig, end1, neighbor, contig.otherEndOfLinks[end1][n])
                                        
                            # print("RIght of 138 there is ", [i.ID for i in segments[names['138']].links[1]], " ", [i for i in segments[names['119']].otherEndOfLinks[1]], " ", segments[names['119']].ID)
                            # print("RIght of 119 there is ", [i.ID for i in segments[names['119']].links[1]], " ", [i for i in segments[names['138']].otherEndOfLinks[1]], " ", segments[names['138']].ID)
                                
                            if len(contig.links[1-end1]) > 0 and c < len(contigs)-1: #delete all the links right of the contig, the only good one will be reestablished later
                                while len(contig.links[1-end1]) > 0 :
                                    neighbor = contig.links[1-end1][0]
                                    segment.delete_link(contig, 1-end1, neighbor, contig.otherEndOfLinks[1-end1][0])

                            
                            #now the link the contig to the contig right at its left
                            segment.add_link(contig, end1, segments[newContigsIndices[c-1]], end0, CIGAR)                                
                        
                            newContigsIndices += [oldContigsIndices[c]]
                            
                        CIGAR = nextCIGAR

                                                
                    