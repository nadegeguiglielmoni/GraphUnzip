#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:19:55 2020

File dedicated to the analysis of the HiC coverage of contigs
"""

import matplotlib.pyplot as plt
import numpy as np
import basic_functions as bf

def determine_HiC_coverage(hiccontactsfile, info_contig, fragment_list) : #returns the number of HiC contacts per basepair
    
    with open(hiccontactsfile) as f:
    
        coverage = [0]*(info_contig[-1][0]+1)

        for line in f :
            
            line = line.strip('\n')
            line = line.split('\t')
            
            if line[0] != '487796' :#because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]
                
                contig1 = fragment_list[contact[0]][0]
                contig2 = fragment_list[contact[1]][0]
                coverage[contig1] += contact[2]
                coverage[contig2] += contact[2]
    
    xrange = [i for i in range(len(coverage))] 
    plt.scatter(xrange, coverage)
    #plt.ylim([0,1])
    
    return coverage

def determine_unconnected_contigs(hiccontactsfile, fragmentList) :
    
    #first we are going to look at which contigs have no HiC contacts
    with open(hiccontactsfile) as f:
    
        does_a_contig_interact_with_this_one = [False for i in range(fragmentList[-1][0]+1)]
        for line in f :
            
            line = line.strip('\n')
            line = line.split('\t')
            
            if line[0] != '487796' :#because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]
                   
                contig1 = fragmentList[contact[0]][0]
                contig2 = fragmentList[contact[1]][0]
            
                does_a_contig_interact_with_this_one[contig1] = True
                does_a_contig_interact_with_this_one[contig2] = True
                    
        contigUnconnected = []
        for i in range(len(does_a_contig_interact_with_this_one)):
            if does_a_contig_interact_with_this_one[i] == False :
                contigUnconnected += [i]
                
        print(contigUnconnected)
 
def check_if_there_are_restriction_fragments_in_this_contig(contig, restrictionSiteSequence, genomeFastaFile) :
    
    with open(genomeFastaFile) as f :
        
        target = '>sequence'+str(contig)
        readnextline = False
        
        for line in f :
            
            if readnextline == True :
                                
                return line.count(restrictionSiteSequence)
            else :
                #print(target)
                if target in line :
                    readnextline = True
      
    print('There is a problem with the input contig')                              
    return 0

def check_if_there_are_restriction_fragments_in_unconnected_contigs(contigs, restrictionSiteSequence, genomeFastaFile) :
    
    for i in contigs :
        print ('The contig ' + str(i) + ' contains '+ \
               str(check_if_there_are_restriction_fragments_in_this_contig(i, restrictionSiteSequence, genomeFastaFile))+\
                   ' restriction sites')

def correlation_GCcontent_HiCcoverage(coverage, genomeFastaFile, unconnectedContigs):
    
    GCcontentOfContigs = [-1]*len(coverage)
    
    with open(genomeFastaFile) as f :
    
        step = 0
        for line in f :
            
            if not '>' in line :
                GCcontent = (line.count('G')+line.count('C')) /len(line)
                GCcontentOfContigs[int(step/2)] = GCcontent
                    
            step += 1
        
        extract = [GCcontentOfContigs[i] for i in range(len(GCcontentOfContigs)) if i in unconnectedContigs]
        extractCoverage = [coverage[i] for i in range(len(coverage)) if i in unconnectedContigs]
        #plt.scatter(extract, extractCoverage)
        plt.scatter(GCcontentOfContigs, coverage)
        plt.xlabel('GC content')
        plt.ylabel('Coverage')
        plt.ylim([0,1])
    

fragmentList = bf.import_from_csv('listsPython/fragmentList.csv')
infcontigs = bf.read_info_contig('data/results/info_contigs.txt')
#unconnectedcontigs = bf.import_from_csv('listsPython/unconnectedContigs.csv')
#unconnectedcontigs = [x[0] for x in unconnectedcontigs]
#print(len(unconnectedcontigs))
#check_if_there_are_restriction_fragments_in_unconnected_contigs(unconnectedcontigs, 'GATC','data/Assembly.fasta') #GATC corresponds to the cutting site of DpnII

coverage = determine_HiC_coverage('data/results/abs_fragments_contacts_weighted.txt', infcontigs, fragmentList)
bf.export_to_csv(coverage, 'listsPython/HiCcoverage.csv')
#coverage = bf.import_from_csv('listsPython/HiCcoverage.csv')
#coverage = [x[0] for x in coverage]
    
#coverageWithoutTheZeros = [x for x in coverage if x != 0]
#coverageLog = [np.log10(x) for x in coverageWithoutTheZeros]

#plt.hist(coverageLog, bins = 160)
#plt.xlabel('Number of HiC contacts per bp (log scale)')
#plt.ylabel('Number of contigs')
#plt.xlim([0,3])

#correlation_GCcontent_HiCcoverage(coverage, 'data/Assembly.fasta', unconnectedcontigs)
print('Finished')