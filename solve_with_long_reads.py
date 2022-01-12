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
from finish_untangling import merge_adjacent_contigs
from input_output import read_GAF
from input_output import read_TSV
import time

import segment

#Master function of the file
#Input : initial gfa (as a list of segments), a GAF file with long reads mapped to the segments, names (which is an index numbering the contig), multiplicity : the pre-computed ploidy of each contig (as numbered in names)
#Output : new gfa (as a list of segments) corrected with long reads, and modified copiesnumber (taking into account contigs that have been duplicated)
def bridge_with_long_reads(segments, names, copiesnumber, gafFile, supported_links2, refHaploidy, multiplicities):
        
    ##There are two phases : first build consensus with approximate haploid contigs. Detect the incoherence, obtain a reliable haploid list and do that all over again
    
    #first phase : determine all contigs that look haploid, and sort them by length
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file    
    #print("multiplicity of 95 : ", multiplicities[names['edge_95']])
    #print(multiplicities)
    # haploidContigs = []
    # for se, s in enumerate(segments) :
        
    #     if multiplicities[se] == 1 :
    #     #to be deemed haploid, a segment must have at most one connection at each of its end plus be equally or less covered thant neighboring contigs
    #     #if depth could not be read, all contigs will be supposed haploid and we'll iteratively correct this assumption in a later phase
    #         links = s.links
    #         if round(s.depths[0]/refHaploidy) == 1 :
    #             haploidContigs.append(s)
    #         elif len(links[0]) == 1 and len(links[1]) == 1 :# and s.depths[0] < 1.5*links[0][0].depths[0] and s.depths[0] < 1.5*links[1][0].depths[0] : 
    #             haploidContigs.append(s)
    #         elif len(links[0]) == 1 and len(links[1]) == 0  : #and s.depths[0] < 1.5*links[0][0].depths[0] : 
    #             haploidContigs.append(s)
    #         elif len(links[0]) == 0 and len(links[1]) == 1 : # and s.depths[0] < 1.5*links[1][0].depths[0] : 
    #             haploidContigs.append(s)
        
    
    # haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    # haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    # index = 0
    # for s in haploidContigs :
    #     haploidContigsNames[s.names[0]] = index
    #     index += 1
    
    
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file
    
    lines = []
    if '.gaf' in gafFile :
        print("Reading the gaf file...")
        read_GAF(gafFile, 0.7, 0.1, lines)
        print("Finished going through the gaf file.")
    elif '.tsv' in gafFile :
        print("Reading the tsv file...")
        read_TSV(gafFile, names, lines)
        print("Finished going through the tsv file.")
    else :
        print("ERROR: input format of mapped read not recognized. It should be .gfa or .gpa")
        sys.exit()
    
    #determine an approximate list of contigs that look haploid
    haploidContigs, haploidContigsNames = determine_haploid_contigs(lines, segments, names)
    #print(haploidContigsNames)
    sure_haploids = False
    
    #inventoriate all bridges in the list bridges : sequence of contigs found in the GAF containing at least one contig of multiplicity 1. Do fill in the longContigs list while you're at it
    longContigs = [] #long contigs are contigs too long to be traversed by a long read in the gaf
    longContigs = [True for i in range(len(names))] #then all contigs that are in the middle of a read will be marked as False
    bridges = [[[],[]] for i in range(len(haploidContigs))] #bridges is a list inventoring at index haploidCOntigsNames[seg.names[0]] all the links left and right of the contig, supported by the gaf
    minimum_supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #minimum_supported links is the list of all links between different contigs found at least once in the gaf file
    inventoriate_bridges(lines, bridges, minimum_supported_links, haploidContigsNames, longContigs, names)
    
    #now, from all the bridges, build consensus bridges
    consensus_bridges = [['',''] for i in range(len(haploidContigs))] #consensus bridge is essentially the same as bridges, except there is only one bridge left at each side for each contig
    print("Building consensus bridges from all the long reads")

    haploidContigs, haploidContigsNames, consensus_bridges = build_consensus_bridges(consensus_bridges, bridges, names, haploidContigs, haploidContigsNames)
    print("Done building consensus bridges")

    #print("Consensus bridge of contig 223445 : ", consensus_bridges[haploidContigsNames['223445']])
    
    bridges = []
        
    print("Now we will determine through an iterative process what contigs of the assembly are present only once in the final genome")
    while not sure_haploids : #knowing what contigs are really haploid may take several iterations
        
        leng = len(haploidContigs)
        
       # print("consensus of 20: ", consensus_bridges[haploidContigsNames['20']])


        #consensus bridges overlap two by two (e.g. >2>3>4 right of 1 overlaps with <3<2<1 left of 4), so merge them to have a set of non-overlapping consensus bridges
        non_overlapping_bridges = [['',''] for i in range(len(haploidContigs))] 

        haploidContigs, haploidContigsNames, consensus_bridges, sure_haploids = merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs, longContigs, names, multiplicities)

        #print("nonoverlap consensus of 20: ", non_overlapping_bridges[haploidContigsNames['20']])
        #if this last phase of merge_contig detected no inconstistencies, sure_haploids=True and the program moves on
        
        if not sure_haploids :
            print("Out of ", leng, " supposed single-copy contigs, ", leng-len(haploidContigs), " were not actually haploid. Recomputing until all the single-copy contigs are robust")

    #from the consensus bridges, mark all links that are supported by the long reads
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file, and in how many different consensus
    compute_supported_links(supported_links, consensus_bridges, haploidContigsNames, haploidContigs, longContigs, names)
    
    supported_links = supported_links.tocoo()
    supported_links2 = supported_links2.todok()
    #now the haploid contigs are determined with confidence
    
    for r, c, m in zip(supported_links.row, supported_links.col, supported_links.data):
        supported_links2 [r,c] = max(supported_links2[r,c], m)  
        supported_links2 [c,r] = max(supported_links2[c,r], m)  
    #print("Link is supported after maxing with strength ", supported_links2[2*names['127'], 2*names['112']])
    # print("Link is supported after maxing with strength ", supported_links2[2*names['2505'], 2*names['2489']])
         
    #print("Contig 1621 has for consensus : ", consensus_bridges[haploidContigsNames["1621"]])

    #now actually unzip the graph using the instructions in non_overlapping_bridges
    print("Let's move on to actually untangling the graph")
    unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, haploidContigs, haploidContigsNames, names, supported_links2, minimum_supported_links, multiplicities, longContigs)
    
    #now remove the tips that apparently came from an overestimation in the multiplicity
    print("Now we correct the last quirks by looking a posteriori at the graph")
    merge_adjacent_contigs(segments)
    trim_tips(segments, multiplicities, names, haploidContigsNames)
        
    #print(non_overlapping_bridges)
    
    return segments


#input : all aligned long reads
#output : A list of "haploid" contigs, i.e. contigs that have at most one possible other contig right and left
def determine_haploid_contigs(lines, segments, names) :
    
    #haploidContigsIdx = set([i for i in range(len(segments))]) #list of the idx of all haploid contigs : we'll whittle it down
    neighborLeftRight = [({}, {}) for i in range(len(segments))] #list of neighbor contigs left and right of each contig : if a contig has only one neighbor left and right we'll say it's haploid
    
    for line in lines :
        
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''
    
        for c, contig in enumerate(contigs) :
            
            #if names[contig] in haploidContigsIdx :
                
                orientation = '><'.index(orientations[c])
                
                if c > 0 :
                    
                    if contigs[c-1] not in neighborLeftRight[names[contig]][orientation] :
                        neighborLeftRight[names[contig]][orientation][contigs[c-1]] = 1
                    else :
                        neighborLeftRight[names[contig]][orientation][contigs[c-1]] += 1
                    # if len (neighborLeftRight[names[contig]][orientation].keys()) > 1 :
                    #     #haploidContigsIdx.discard(names[contig])
      
                if c < len(contigs) -1 :
                    
                    if contigs[c+1] not in neighborLeftRight[names[contig]][1-orientation] :
                        neighborLeftRight[names[contig]][1-orientation][contigs[c+1]] = 1
                    else :
                        neighborLeftRight[names[contig]][1-orientation][contigs[c+1]] += 1
                        
                    # if len (neighborLeftRight[names[contig]][1-orientation].keys()) > 1 :
                    #     #haploidContigsIdx.discard(names[contig])

    
    haploidContigs = []
    for se, nei in enumerate(neighborLeftRight) :
       
        if (len(nei[0]) ==0 or max(nei[0].values()) > 0.9 * sum(nei[0].values())) and (len(nei[1])==0 or max(nei[1].values()) > 0.9 * sum(nei[1].values())) :
            
            haploidContigs += [segments[se]]
        
       
    haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in haploidContigs :
        haploidContigsNames[s.names[0]] = index
        index += 1
        
    return haploidContigs, haploidContigsNames
        
#input : a list of alignments of a gaf file
#output : the completed bridges list, with for each haploid contig a list of what was found left and right of the contig. 
def inventoriate_bridges(lines, bridges, minimum_supported_links, haploidContigsNames, longContigs, names) :
    
    
    for l, line in enumerate(lines) :      
                    
        if (l+1) % 1000 == 0 :
            print("Inventoried ", l+1, " long reads over ", len(lines))
        
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''
    
        for c, contig in enumerate(contigs) :
            
            if c>0 :
                minimum_supported_links[2*names[contigs[c-1]] + '<>'.index(orientations[c-1]) , 2*names[contigs[c]] + '><'.index(orientations[c])] = 1
                minimum_supported_links[2*names[contigs[c]] + '><'.index(orientations[c]), 2*names[contigs[c-1]] + '<>'.index(orientations[c-1])] = 1
    
            if c > 0 and c < len(contigs) - 1 :
                longContigs[names[contig]] = False
                                        
            if contig in haploidContigsNames :
                
                
                if orientations[c] == ">" :
                    r = 0
                    #first look at what contigs are left of the contig of interest
                    bridges[haploidContigsNames[contig]][1] +=  [""]
                    for c2 in range(c+1, len(contigs)) :
                        
                        bridges[haploidContigsNames[contig]][1][-1] += orientations[c2] + contigs[c2]
                        
                        # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                        #     break
                        
                    #then look at what's left of the contig of interest (so mirror the orientations)
                    bridges[haploidContigsNames[contig]][0] +=  [""]
                    for c2 in range(c-1 , -1, -1) :
                        
                        if orientations[c2] == '>' :
                            bridges[haploidContigsNames[contig]][0][-1] += '<' + contigs[c2]
                        else :
                            bridges[haploidContigsNames[contig]][0][-1] += '>' + contigs[c2]
                            
                        # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                        #     break
                            
                else :
                    
                    #first look at what contigs are left of the contig of interest
                    bridges[haploidContigsNames[contig]][0] +=  [""]
                    for c2 in range(c+1, len(contigs)) :
                        bridges[haploidContigsNames[contig]][0][-1] += orientations[c2] + contigs[c2]
                        # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                        #     break
                        
                    #then look at what's left of the contig of interest (so mirror the orientations)
                    bridges[haploidContigsNames[contig]][1] +=  [""]
                    for c2 in range(c-1 , -1, -1 ) :
                        
                        if orientations[c2] == '>' :
                            bridges[haploidContigsNames[contig]][1][-1] += '<' + contigs[c2]
                        else :
                            bridges[haploidContigsNames[contig]][1][-1] += '>' + contigs[c2]
                            
                        # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                        #     break
                                                 
#input : list of bridges for each haploid contig
#output : completed consensus_bridges, where there is max one bridge at each end of contig
def build_consensus_bridges(consensus_bridges, bridges, names, haploidContigs, haploidContigsNames):
    
    not_actually_haploid = [] #a list of not actually haploid contigs : at this point, we will rule out as haploid contigs those that consensus back on themselves
    
    #print(allcontigs[0][1])
    
    for c in range(len(bridges)) :
                
        if (c)%100 == 0 :
            print("consensused ", c, " bridges out of ", len(consensus_bridges), end='\r')
            
        localContigs = [ [ re.split('[><]' , bridges[c][j][k])[1:] for k in range(len(bridges[c][j])) ] for j in range(2)]
        localOrientations = [ [ "".join(re.findall("[<>]", bridges[c][j][k])) for k in range(len(bridges[c][j])) ] for j in range(2)]
        
        for end in range(2) : #look left and right of each contig
        
            kept_reads = [i for i in range(len(bridges[c][end])) if bridges[c][end][i] != ""] #list keeping in mind all the reads still in the consensus
            pos = 0
            consensus = ''

            while len(kept_reads) > 0 :
                
                candidate2 = [localContigs[end][i][pos] for i in kept_reads]
                candidate1 = [localOrientations[end][i][pos] for i in kept_reads]
                                
                cons1 = Counter(candidate1).most_common(1)[0]
                cons2 = Counter(candidate2).most_common(1)[0]
                                
                if (cons1[1] > 0.8*len(kept_reads) or (cons1[1]>=2 and cons1[1] == len(kept_reads)-1)) and (cons2[1] > 0.8*len(kept_reads) or (cons2[1]>=2 and cons2[1] == len(kept_reads)-1)) : #consensus is defined there as a 80% of the propositions agreeing or only 1 proposition disagreeing
                
                    if cons1[0] == '*' :
                        break
                    
                    consensus_bridges[c][end] = consensus_bridges[c][end] + cons1[0] + cons2[0]
                    
                    # if cons2[0] == haploidContigs[c].names[0] and len(consensus_bridges[c][end]) < len(names): #if the contig loops on itself, it is not haploid except if we're on a circular chromosome
                    #     #not_actually_haploid += [c]
                    #     #break
                    #     pass
                    
                    
                    #inventoriate this bridge in supported_links
                    
                    #print('Adding from ', consensus_bridges[c][end], " c: ", c)
                    
                    #update the list of reads still there to build the consensu
                    new_kept = []
                    for r in kept_reads:
                        if len(localContigs[end][r]) > pos+1 :
                            if localOrientations[end][r][pos] == cons1[0] and localContigs[end][r][pos] == cons2[0] :
                                new_kept += [r]
                        
                    kept_reads = new_kept
                    
                    # if cons2[0] in haploidContigsNames : #just go up to the next haploid contig
                    #     break
               
                else : #if there is no consensus let's stop
                    # if "11832294" in haploidContigs[c].names :
                    #     print("There are no consensus there: ", candidate1, candidate2)
                    break
                
                pos += 1
        bridges[c] = [] #to free memory
        
    new_consensus_bridges = []
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
            new_consensus_bridges += [consensus_bridges[i]]
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
                    
    return reliable_haploid_contigs, reliable_haploid_contigsNames, new_consensus_bridges
                
#input: list of consensus bridges
#output: filled supported_links matrix
def compute_supported_links(supported_links, consensus_bridges, haploidContigsNames,haploidContigs, longContigs, names) :
    
    for c in range(len(consensus_bridges)) :
        
        for end in range(2) :
            
            
            contigs = re.split('[><]' , consensus_bridges[c][end])
            orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
            del contigs[0] #because the first element is always ''
            
            previousContig = names[haploidContigs[c].names[0]] # variable containing the index of the last contig of the consensus (useful to fill supported_links)
            previousEnd = end #variable containing the free end of the last contig of the consensus (useful to fill supported_links)
            firstHapContigMet = False #while no haploid contig is in the consensus, links should be supported (that's the way to ensure all supported links are supported only once)
            
            nobridges = True #to see if there is another end to the bridge
            for contig in contigs:
                if contig in haploidContigsNames :
                    nobridges = False
            
            halfbridge = False
            if len(contigs) > 0 :
                halfbridge = (nobridges and longContigs[names[contigs[-1]]]) #to mark where there are half bridges, where links are supported both ways
            
            for co, contig in enumerate(contigs):

                # if c == haploidContigsNames['83467'] or c == haploidContigsNames['55374'] :

                current_end = 0
                if orientations[co] == "<" :
                    current_end=1
                if not firstHapContigMet :
                    supported_links[2*names[contig]+current_end, 2*previousContig+previousEnd] += 1  #each link is present twice in supported_links, only add to one cell of supported_links, the symetric will be added in bridge_with_long_reads
                    if halfbridge :
                        supported_links[2*previousContig+previousEnd, 2*names[contig]+current_end] += 1 
                      
                previousContig = names[contig]
                previousEnd = 1-current_end
                
                if contig in haploidContigsNames and not firstHapContigMet:
                    firstHapContigMet = True

                    # if c == haploidContigsNames['83467'] or c == haploidContigsNames['55374'] :
                    #     print("First hap contig encountered: ", contig)
                                                        
#input: list of all consensus bridges for all haploid contigs, plus a list of the haploid contigs
#output: a list of non-overlapping bridges left and right of each contig, with only full bridges (connecting 2 haploid contigs)
def merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs, longContigs, names, multiplicities) :
    
    sure_haploids = True
    not_actually_haploid = [] #a list of not actually haploid contigs among the haploidContigs
    
    for c in range(len(consensus_bridges)) :
        
        for end in range(2) :
            
            if consensus_bridges[c][end] != '' :
            
                contigs = re.split('[><]' , consensus_bridges[c][end])
                orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
                del contigs[0] #because the first element is always ''
                
                firstHapIdx = -10000
                for co in range(len(contigs)):
                    if contigs[co] in haploidContigsNames :
                        firstHapIdx = co
                        break
                
                if firstHapIdx != -10000 : #in case there is a full bridge
                    #check if the two bridges are coherent
                    coherent = False
                    otherEnd = 0
                    if orientations[firstHapIdx] == "<":
                        otherEnd = 1
                        
                        # print ("here is contigs [0] : ", haploidContigs[c].names[0], " -- ", re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[-1]]][otherEnd]))
                    
                    symmetrical = re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[firstHapIdx]]][otherEnd])
                    
                    if len(symmetrical)>firstHapIdx+1 and haploidContigs[c].names[0] == symmetrical[firstHapIdx+1] : #+1 because the first element of symmetrical is always ''
                        coherent = True
                        
                    if not coherent :
                        
                        #check if the incoherence is not due to being on a very weakly supported (and probably false) path
                        # if all([i not in haploidContigsNames for i in symmetrical]) or True: #else, it means that the path is very minoritary compared to the path going to the contig making this assumption false
                        
                            #print("two incoherent bridges between ", haploidContigs[c].names[0] , " and ", haploidContigs[haploidContigsNames[contigs[firstHapIdx]]].names[0], " : ", consensus_bridges[c][end], " vs ", consensus_bridges[haploidContigsNames[contigs[firstHapIdx]]][otherEnd])
                        not_actually_haploid += [haploidContigsNames[contigs[firstHapIdx]]]
                        sure_haploids = False
                            
                            # print("Contig ", contigs[firstHapIdx], " does not look haploid, seen from contig ",  haploidContigs[c].names[0], " : ", contigs[:firstHapIdx+1], ", whose bridge is ", symmetrical[:firstHapIdx+2])

                        # else : #if the path is really minoritary, the contig is probably not haploid, or at least you should not count on it
                        #     not_actually_haploid += [c]
                        #     sure_haploids = False
                        #     if haploidContigs[c].names[0] == '7' :
                        #         print("Contig ", haploidContigs[c].names[0], " does not look haploid, here is its bridge to the next contig ", contigs[firstHapIdx] , " : ", contigs[:firstHapIdx+1], ", whose bridge is ", symmetrical[:firstHapIdx+2])
                        #     #print("Haploid contigs : ", haploidContigsNames)
                            
                        
                        
                    else :
                        if c < haploidContigsNames[contigs[firstHapIdx]] : #keep only one of the two bridges
                            non_overlapping_bridges[c][end] = ''.join([orientations[i]+contigs[i] for i in range(firstHapIdx+1)])
                
                elif longContigs[names[contigs[-1]]] : #if it ends on a  multiploid long contig, extracting could still be done (though not for the last, long contig). That's because we know there is no bridge going the other direction through the long contig
                    non_overlapping_bridges[c][end] = consensus_bridges[c][end]
                
                else : #a half-bridge cannot be unambiguously extacted from the graph
                    non_overlapping_bridges[c][end] = ''
    
    #now update haploidCOntig and consensus_bridges
        
    new_consensus_bridges = []
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    not_actually_haploid = list(set(not_actually_haploid))
    not_actually_haploid.sort()
        
    index = 0
    indexNot = 0
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(not_actually_haploid) and i == not_actually_haploid[indexNot] :
            indexNot += 1
            if multiplicities[names[haploidContigs[i].names[0]]] == 1 :
                multiplicities[names[haploidContigs[i].names[0]]] += 1
        else :
            new_consensus_bridges += [consensus_bridges[i]]
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
                    
    return reliable_haploid_contigs, reliable_haploid_contigsNames, new_consensus_bridges, sure_haploids #return the updated list of haploid contigs
          
                
#input : a list of segments and the non_overlapping_bridges
#output: a list of segment where all the bridges have been built          
def unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, haploidContigs, haploidContigsNames, names, supported_links, minimum_supported_links, multiplicities, longContigs) :  

    #compute the minimum multiplicity of each contig, so that there is enough contig to build all bridges, even when the depth suggests otherwise
    #first check with supported links
    for s, seg in enumerate(segments) :
        minLeft = 0
        minRight = 0
        for n, neighbor in enumerate(seg.links[0]) :
            minLeft += supported_links[2*names[seg.names[0]], 2*names[neighbor.names[0]]+seg.otherEndOfLinks[0][n]]
        for n, neighbor in enumerate(seg.links[1]) :
            minRight += supported_links[2*names[seg.names[0]]+1 , 2*names[neighbor.names[0]]+seg.otherEndOfLinks[1][n]]
        
        multiplicities[s] = max( min (minLeft, minRight), multiplicities[s])
    #second check with non_overlapping_bridges
    minimum_multiplicity = [0 for i in range(len(multiplicities))] 
    for c in range(len(non_overlapping_bridges)) :
        for end in range(2) :
            contigs = re.split('[><]' , non_overlapping_bridges[c][end])
            del contigs[0] #because the first element is always ''
            contigs = contigs[:-1]
            for contig in contigs :
                minimum_multiplicity[names[contig]] += 1
    for s in range(len(segments)) :
        multiplicities[s] = max(minimum_multiplicity[s], multiplicities[s])
            
    alreadyDuplicated = [-1 for i in range(len(names))] #list useful for duplicating long contigs : only duplicate them from one side, the one that is in this list 
            
    l = len(segments)
    for se in range(l) :
        
        if (se)%1000 == 0 :
            print("Processed ", se, " contigs out of ", l, ", while untangling with long reads")
        s = segments[se]
         
        if s.names[0] in haploidContigsNames :
            
                        
            for end in range(2) :
                
                
                if len(non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]) > 0  :# and ('2' not in s.names or end == 1) and ('7' not in s.names or end == 0 ): #means there is a bridge to be built there
                
                    #print("Unzipping contig ", s.names[0], " with bridge : ")
                                                    
                    contigs = re.split('[><]' , non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                    orientations = "<>"[end] + "".join(re.findall("[<>]", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]))
                    del contigs[0] #because the first element is always ''
                    contigs = [s.names[0]] + contigs
                    contigsID = [segments[names[c]].ID for c in contigs]
                    newContigsIndices = [names[contigs[0]]] #list keeping track of the indices of the contigs on the path
                    oldContigsIndices = [names[i] for i in contigs]
                                        
                    haploidCoverage = (segments[names[contigs[0]]].depths[0]+segments[names[contigs[-1]]].depths[0])/2 #computation of the haploid coverage at this point in the assembly
                    
                    nextEnd = 0
                    if orientations[1] == '<' :
                        nextEnd = 1
                    CIGAR = segments[names[contigs[0]]].CIGARs[end][segment.find_this_link(segments[oldContigsIndices[1]], nextEnd, segments[names[contigs[0]]].links[end], segments[names[contigs[0]]].otherEndOfLinks[end])]
                    nextCIGAR = '';

                    
                    #take care of the first contig, nobody is going to do it elsewhise
                    contig = segments[names[contigs[0]]]
                    idx = 0
                    if len(contig.links[end]) > 0 : #delete all the links right of the contig, the only good one will be reestablished later (actually, only delete the links to haploidContigs)
                        #print('here edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)
                        while len(contig.links[end]) > idx :
                            neighbor = contig.links[end][idx]
                            if neighbor.names[0] in haploidContigsNames or alreadyDuplicated[names[neighbor.names[0]]] != 1-nextEnd :
                                success = segment.delete_link(contig, end, neighbor, contig.otherEndOfLinks[end][idx])
                            else :
                                idx += 1

                        #print('now edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)

                    #then take care of all other contigs
                    for c in range(1, len(contigs)) :
                        
                        contig = segments[names[contigs[c]]]
                            
                        #multiplicity = max([1, round(contig.depths[0]/haploidCoverage), multiplicities[names[contigs[c]]]]) #the multiplicity is inferred from the coverage of the contig compared to the coverage of the two haploid contigs at the extremity of the bridge
                        multiplicity = multiplicities[names[contigs[c]]]
                        multiplicities[names[contigs[c]]] -= 1
                        # if contigs[c] == '2505' :
                        #     print("Using once contig 2505 in overlap ", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end], ", now at : ", multiplicities[names[contigs[c]]])
                        #neighborMultiplicity = max([1, round(segments[oldContigsIndices[c-1]].depths[0]/haploidCoverage), multiplicities[names[contigs[c-1]]]])
                        neighborMultiplicity = multiplicities[names[contigs[c-1]]]

                        #print(alreadyDuplicated[names['edge_291']])
                        
                        end1, end0 = 1, 1
                        if orientations[c] == '>' :
                            end1 = 0
                        if orientations[c-1] == '<' :
                            end0 = 0
                            
                        #remember the CIGAR before deleting all the links
                        if c < len(contigs)-1 :
                            nextEnd = 0
                            if orientations[c+1] == '<' :
                                nextEnd = 1
                            if segment.find_this_link(segments[oldContigsIndices[c+1]], nextEnd, contig.links[1-end1], contig.otherEndOfLinks[1-end1]) != -1 :
                                nextCIGAR = contig.CIGARs[1-end1][segment.find_this_link(segments[oldContigsIndices[c+1]], nextEnd, contig.links[1-end1], contig.otherEndOfLinks[1-end1])]
                            else :
                                print("Debug WARNING, ", contigs, " : looking for ", segments[oldContigsIndices[c+1]].names, " ", nextEnd, " from ", contig.names, " among ", [i.names for i in contig.links[1-end1]], " ", contig.otherEndOfLinks[1-end1], " ", s.names[0], " ",non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                            
                        #print("Link supported with strength ", supported_links[names[contigs[c]*2+end1 , names[contigs[c-1]]*2+end0])
                        
                        # if s.names[0] == '55374' and end == 0 and contigs[c] == '2505' :
                        #     print("multiplicity of 2505 : " , multiplicity)
                                                    
                        if multiplicity > 1 and (c < len(contigs)-1 or (longContigs[names[contigs[-1]]] and alreadyDuplicated[names[contigs[-1]]] != 1 - end1)): #if multiplicity>1, the contig should be duplicated 
                        
                            newSegment = segment.Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i/multiplicity for i in contig.depths])
                            segments.append(newSegment)
                            newContigsIndices += [len(segments) - 1]
                            
                            for n in newSegment.names :
                                copiesnumber[n] += 1

                            #add the link to form the new bridge
                            segment.add_link(segments[-1], end1, segments[newContigsIndices[c-1]], end0, CIGAR)

                            #delete the old link if and only if it was only supported by one path only
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            supported_links[names[contigs[c-1]]*2+end0, names[contigs[c]]*2+end1] -= 1
                            minimum_supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            minimum_supported_links[names[contigs[c-1]]*2+end0, names[contigs[c]]*2+end1] -= 1
                                
                            if supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] == 0 and len(segments[oldContigsIndices[c-1]].links[end0])>1 and segment.find_this_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]].links[end0], segments[oldContigsIndices[c-1]].otherEndOfLinks[end0], warning=False) != -1:
                                # if "130" in segments[oldContigsIndices[c-1]].names : 
                                #     print("remove links left of ", segments[oldContigsIndices[c]].names, " : ", segments[oldContigsIndices[c-1]].names)
                                segment.delete_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]], end0, warning = True) #though it is very hard to be sure, it is not impossible at that point that we actually delete a link that is present in another bridge, so don't warn
                                
                            #since contig has been duplicated, lower its depth
                            contig.divide_depths(multiplicity/(multiplicity-1))
                            
                            if c == len(contigs) - 1 :
                                for n, neighbor in enumerate(contig.links[1-end1]) :
                                    segment.add_link(segments[-1], 1-end1, neighbor, contig.otherEndOfLinks[1-end1][n], contig.CIGARs[1-end1][n])
                                alreadyDuplicated [names[contigs[-1]]] = end1
                                    
                                    
                        elif c == len(contigs)-1 and longContigs[names[contigs[c]]] and alreadyDuplicated[names[contigs[-1]]] == 1 - end1 : #if the contig has already been duplicated from the other side 

                            if newContigsIndices[-1] != oldContigsIndices[c-1]: #only need to add links if you just duplicated

                                for alreadyDuplicatedContig in segments[oldContigsIndices[c-1]].links[end0] :
                                    if alreadyDuplicatedContig.names[0] == contigs[c]:
                                        segment.add_link(segments[newContigsIndices[-1]] , end0, alreadyDuplicatedContig, end1, CIGAR)

                                                        
                        else :
                            

                            # if the contig has multiple incoming reads left, get rid of all of them (the good one is reestablished later). Do the same at the other end of the contig
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            minimum_supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1       
                            
                            if len(contig.links[end1]) > 0 :
                                for n, neighbor in enumerate(contig.links[end1]) :
                                    if neighbor.ID != segments[newContigsIndices[c-1]].ID:
                                        if supported_links[names[contigs[c]]*2+end1 , names[neighbor.names[0]]*2+contig.otherEndOfLinks[end1][n]] <= 0 and minimum_supported_links[names[contigs[c]]*2+end1 , names[neighbor.names[0]]*2+contig.otherEndOfLinks[end1][n]] <= 0  :
                                            segment.delete_link(contig, end1, neighbor, contig.otherEndOfLinks[end1][n])
                                        
                            # print("RIght of 138 there is ", [i.ID for i in segments[names['138']].links[1]], " ", [i for i in segments[names['119']].otherEndOfLinks[1]], " ", segments[names['119']].ID)
                            # print("RIght of 119 there is ", [i.ID for i in segments[names['119']].links[1]], " ", [i for i in segments[names['138']].otherEndOfLinks[1]], " ", segments[names['138']].ID)
                            if len(contig.links[1-end1]) > 0 and c < len(contigs)-1 : 
                                
                                nextEnd = 2
                                for n, neighbor in enumerate(contig.links[1-end1]) :
                                    if neighbor.names[0] == contigs[c+1] :
                                        nextEnd = contig.otherEndOfLinks[1-end1][n]
                                
                                if not (c == len(contigs)-2 and alreadyDuplicated[names[contigs[c+1]]] == 1-nextEnd and contigs[c+1] not in haploidContigs): #delete all the links right of the contig, the only good one will be reestablished later, except if looking at a long contig duplicated from the other end
                                    while len(contig.links[1-end1]) > 0 :
                                        neighbor = contig.links[1-end1][0]
                                        segment.delete_link(contig, 1-end1, neighbor, contig.otherEndOfLinks[1-end1][0])

                            
                            #now the link the contig to the contig right at its left
                            segment.add_link(contig, end1, segments[newContigsIndices[c-1]], end0, CIGAR) 
                            
                            newContigsIndices += [oldContigsIndices[c]]
                            
                        CIGAR = nextCIGAR
                    
                        
                    #print('edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)
                # if supported_links[names['2484']*2+1 , names['2505']*2] == 1:
                #     print(non_overlapping_bridges[haploidContigsNames[s.names[0]]][end], "  ", [i.names for i in segments[names['2484']].links[1]], " ", s.names[0], " ", end)
                #     while True :
                #         rien = 0

      
#input : the graph of segments as well as the list of multiplicities
#output : the graph where the dubious tips have been deleted (those where multiplicity has been overestimated)
def trim_tips(segments, multiplicities, names, haploidContigsNames):
    
    toDelete = []
    for s, seg in enumerate(segments):
        
        for end in range(2) :
            
            if len(seg.links[1-end]) == 0 and len(seg.links[end]) == 1 :

                #print("Checking if ", seg.names, " is a tip")
                if any([extended_length(i, seg.links[end][0].otherEndOfLinks[seg.otherEndOfLinks[end][0]][e], 10*seg.length, 30) for e,i in enumerate(seg.links[end][0].links[seg.otherEndOfLinks[end][0]])]) : #this means we're in a very short dead end
                                    
                    if all([i not in haploidContigsNames for i in seg.names]) : #then it means it's probably an error in the determination of the multiplicity
                    
                        segment.delete_link(seg, end, seg.links[end][0], seg.otherEndOfLinks[end][0])
                        toDelete += [s]
    
    for i in toDelete[::-1]:
        del segments[i]

#a function returning True if you can go far (up to threshold) with neighbors of neighbors of neighbors of... 
#it returns False if it needs to recur deeper than thresholdContigs (even though it might be true)
def extended_length(segment, end, thresholdLength, thresholdContigs) :
    
    #print("Extended length called with threshold ", thresholdLength, " on segment , ", segment.names)
    
    if thresholdContigs == 0 :
        return False
    
    if segment.length > thresholdLength :
        return True
    
    #start by looking down the longest contig, it will be fastest
    longestContig = [i for i in range(len(segment.links[1-end]))]
    longestContig.sort(key= lambda x : segment.links[1-end][x].length, reverse = True)
    
    #print("Longest contigs : ", longestContig, [segment.links[1-end][i].length for i in longestContig])
    
    for n in longestContig[:min(len(longestContig), 2)] : #only explore the 2 most promising neighbors, beyond it's not worth it
        
        neighbor = segment.links[1-end][n]
        if extended_length(neighbor, segment.otherEndOfLinks[1-end][n], thresholdLength-segment.length, thresholdContigs-1) :
            return True
        
    return False
        
        
        
        
        
        
        
        
        
        
        
        
        