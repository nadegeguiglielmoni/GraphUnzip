#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 13:02:36 2021

@author: rfaure
"""

#main function of the file : tries to estimate how many copies of each contig is actually present in the actual assembly
#input : a gfa (as a list of segments), with mandatory coverage information ; names (to know at what index to put each contig)
#output : computed_multiplicity, which is a list containing the theoretical multiplicity of contig 'a' at position names[a]
def determine_multiplicity(segments, names) :

    computed_multiplicity = [0 for i in range(len(segments))]
    
    #first compute the 'reference haploid coverage', the average coverage for all reads with axactly one neighbors at each end
    refCoverages = 0
    weightedNumberOfRefContigs = 0
    for s in segments :
        
        links = s.get_links()
        if len(links[0]) == 1 and len(links[1]) == 1 : #if the contig has strictly 1 neighbor at each end we can suppose it is haploid
        
            weightedNumberOfRefContigs += s.length
            refCoverages += s.length * s.depths[0]

      
    refCoverage = refCoverages / weightedNumberOfRefContigs
    
    #then inventoriate all haploid contigs
    for s in segments :
        
        links = s.get_links()
        if len(links[0]) <= 1 and len(links[1]) <= 1 and round(s.depths[0] / refCoverage) == 1 : #if the contig has strictly 1 neighbor at each end we can suppose it is haploid
        
            computed_multiplicity[names[s.names[0]]] = 1
            
    #print("multiplicity of 188 : ", computed_multiplicity[names['edge_188']])
            
       
    #now infer greedily the coverage of contigs next to haploid ones

    while any([computed_multiplicity[i]== 0 for i in range(len(computed_multiplicity))]) :
        unchanged = 0 #keeps the number of unchanged multiplicity since last one
        s = 0 #s is the index of the segment we're looking at
        while unchanged < len(segments) :
            
            if computed_multiplicity[s] == 0 :
                
                new_multiplicity = 0
                new_multiplicity1 = 0
                new_multiplicity2 = 0
                confidence = False
                if len(segments[s].links[0]) > 0 and all([computed_multiplicity[names[i.names[0]]] > 0 for i in segments[s].links[0]]) and all([len(i.links[segments[s].otherEndOfLinks[0][j]]) ==1 for j,i in enumerate(segments[s].links[0])]) :
                    new_multiplicity1 = sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[0]])
                            
                if len(segments[s].links[1]) > 0 and all([computed_multiplicity[names[i.names[0]]] > 0 for i in segments[s].links[1]]) and all([len(i.links[segments[s].otherEndOfLinks[1][j]]) ==1 for j,i in enumerate(segments[s].links[1])]) :
                    new_multiplicity2 = sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[1]])
    
                if new_multiplicity1 == new_multiplicity2 : #in this case the confidence is very high
                    new_multiplicity = new_multiplicity1
                    confidence = True
                else :
                    new_multiplicity = max(new_multiplicity1, new_multiplicity2) #that includes the case where both are > 0, in which case we are not confident
                    
                    
                if new_multiplicity > 0 and (segments[s].depths[0]/refCoverage > new_multiplicity/1.5 or confidence):
                        computed_multiplicity[s] = new_multiplicity
                        #print("coucou, ", segments[s].names, " ", new_multiplicity)
                        unchanged = -1
                # elif new_multiplicity > 0 :
                #     print("not believing coverage ", new_multiplicity, " for contig ", segments[s].names)
                    
            else :
                #check if we can infer the multiplicity of a neighborign contig (possible if there is strictly 1 that is unknown)
                for end in range(2) :
                    
                    if all([len(i.links[segments[s].otherEndOfLinks[end][j]]) ==1 for j,i in enumerate(segments[s].links[end])])  :
                        number0 = 0
                        index0 = 0
                        for n, neighbor in enumerate(segments[s].links[end]) :
                            if computed_multiplicity[names[neighbor.names[0]]] == 0 :
                                number0 += 1
                                index0 = n
                        if number0 == 1 :
                            new_multiplicity = computed_multiplicity[s] - sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[end]])
                            
                            if new_multiplicity > 0 and segments[s].depths[0]/refCoverage > new_multiplicity/1.5 :
                                computed_multiplicity[names[segments[s].links[end][index0].names[0]]] = new_multiplicity
                                #print("salut, ", s+1, " ", names[segments[s].links[end][index0].names[0]]+1, " ", new_multiplicity)
                                unchanged = -1
                            # else :
                            #     print("not believing at all coverage ", new_multiplicity, " for contig ", segments[s].links[end][index0].names[0], " from ", segments[s].names, ". ", computed_multiplicity[s], " and links : ", [computed_multiplicity[names[i.names[0]]] for i in segments[s].links[end]])
                
            unchanged += 1
                    
            
            s = (s+1)%len(segments)
        
        #now the propagation has stopped, start it again by giving a multiplicity to the largest contig without multiplicity yet
        maxi = -1
        maxIdx = -1
        for se, s in enumerate(segments) :
            if computed_multiplicity[se] == 0 :
                if s.length > maxi :
                    maxi = s.length
                    maxIdx = se
        
        if maxi > 0 :

            #print("Adding : ", segments[maxIdx].names, " ", computed_multiplicity[maxIdx], " ",  max(1,round(segments[maxIdx].depths[0] / refCoverage)))
            #computed_multiplicity[maxIdx] = max(1,round(segments[maxIdx].depths[0] / refCoverage))
            minLeft = 0
            minRight = 0
            for n, neighbor in enumerate(segments[maxIdx].links[0]) :
                if len(neighbor.links[segments[maxIdx].otherEndOfLinks[0][n]]) == 1 :
                    minLeft += computed_multiplicity[names[neighbor.names[0]]]
            for n, neighbor in enumerate(segments[maxIdx].links[1]) :
                if len(neighbor.links[segments[maxIdx].otherEndOfLinks[1][n]]) == 1 :
                    minRight += computed_multiplicity[names[neighbor.names[0]]]
                
            computed_multiplicity[maxIdx] = max(1, minLeft, minRight) #choosing 1 there means that computed_multiplicity is a minimum multiplicity

    return refCoverage, computed_multiplicity