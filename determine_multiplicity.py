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
        
        if s.depths[0] == 0 :
            return -1
      
    refCoverage = refCoverages / weightedNumberOfRefContigs
       
    #then infer the multiplicity of each contig by dividing its coverage by the reference haploid covergae
    for s in segments :
        for n, na in enumerate(s.names) : #all sub-contigs have their theoretical multiplicity
            computed_multiplicity[names[na]] = round(s.depths[n] / refCoverage)
    
    # print("The reference coverage is ", refCoverage)
    # for i in range(len(segments)) :
    #     print( segments[i].names[0] , ' has for multiplicity : ' , computed_multiplicity[names[segments[i].names[0]]] )
        
    return refCoverage