#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:05:13 2020

In this file, test functions to test our algorithm
"""

import matplotlib.pyplot as plt
import random
import numpy as np
import time

import basic_functions as bf
from transform_gfa import load_gfa
from basic_functions import export_to_GFA
from solve_ambiguities import solve_ambiguities
from solve_ambiguities import intensity_of_interactions

from evaluate_solution import score_output
from evaluate_solution import draw_distance_HiCcontacts_correlation
from evaluate_solution import heat_solution
from evaluate_solution import simulated_annealing
from loops import flatten_loop

from copy import deepcopy
import scipy.integrate as integrate
from scipy import sparse

def testRatios():
    file = open('ratio.txt', 'r')
    ra = []
    for l in file :
        ra += [float(l)]
        
    plt.hist(ra, bins = 40)
    plt.xlabel('Ratio between the two choices when the programm has to make a choice')
    plt.ylabel('Number of occurences')

def print_chromosomes(chromosomes):
    for c in chromosomes :
        s=''
        for contig in c :
            s += '-'+contig
        print(s[1:])
        
def buildFakeChromosomes(chromosomesLength = 10):
    letters = ['A','A','B','B','C','C','D','D']
    chromosomes = [[letters[j]+str(i) for i in range(chromosomesLength)] for j in range(4)]
    
    duplicates = 0 #number of conitgs that are going to be repeated within each chromosomes :
    for chromosome in chromosomes :
        for i in range(duplicates):
            contigCopied = random.randint(0, len(chromosome)-1)
            insertionPosition = random.randint(0, len(chromosome)-1)
            chromosome.insert(insertionPosition, chromosome[contigCopied])
            
    mutations = 5 #number of total mutations in the genome 
    for i in range(mutations):
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        chromosomes[c][p] += '*'
            
    crossCopies = 1 #number of contigs that are going to be randomly written somewhere in the genome
    for i in range(crossCopies):
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        contigCopies = chromosomes[c][p]
        
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        chromosomes[c].insert(c, contigCopies)
    
    return chromosomes

def exportFakeToGFA(chromosomes, file, lengthOfContig) :
    
    f = open(file, 'w')
    
    alreadyExported = []
    for c in chromosomes :
        for contig in c :
            if contig not in alreadyExported :
                f.write('S\t'+contig+'\t'+ ''.join(['A' for i in range(lengthOfContig)]) +'\n')
                alreadyExported += [contig]
      
    linksAlreadyExported = []
    for c in chromosomes :
        for contig in range(len(c)-1) :
            if [c[contig],c[contig+1]] not in linksAlreadyExported :
                linksAlreadyExported += [[c[contig],c[contig+1]]]
                f.write('L\t'+c[contig]+'\t+\t'+c[contig+1]+'\t+\t100M\n')
                
    return alreadyExported
      
def dist_law(distance) :
    if distance < 10000 :
        distance = 10000
    return 10000/distance

def constructFakeInteractionMatrix(chromosomes, names, lengthOfContigs = 10000):
    
    interactionMatrix = sparse.dok_matrix((len(names), len(names)))
    for c in chromosomes :
        for c1 in range(len(c)):
            for c2 in range(c1, len(c)) :
                con1 = c[c1]
                con2 = c[c2]
                interactionMatrix[names[con1],names[con2]] += integrate.quad(dist_law, (c2-c1-1)*lengthOfContigs, (c2-c1)*lengthOfContigs)[0]
                interactionMatrix[names[con2],names[con1]] += integrate.quad(dist_law, (c2-c1-1)*lengthOfContigs, (c2-c1)*lengthOfContigs)[0]
    
    return interactionMatrix

#function that checks if ls1 is a sublist of ls2 : useful for checking the quality of the output
def sublist(ls1, ls2):
    
    sub = False
    
    s1 = ','.join(ls1)
    s2 = ','.join(ls2)
    sub = sub or (s1 in s2)
    
    s1 = ','.join(ls1[::-1]) #if it goes the other way, it is still imbricated
    s2 = ','.join(ls2)
    sub = sub or (s1 in s2)

    return sub

#function taking as arguments the solution of the problem and the output of the algorithm to see if the output is mistaken
def check_result(chromosomes, listOfSuperContigs, names) :

    #First check if all supercontigs of the output actually exist in the solution (i.e. contigs were not accidentally duplicated)
    for supercontig in listOfSuperContigs :
        
        found = False
        for c in chromosomes :
            if sublist(supercontig.names,c) :
                found = True
        if not found :
            print('Output : ')
            for i in listOfSuperContigs :
                i.print_complete()
            print('Actual chromosomes : ')
            for c in chromosomes :
                print(c)
            print('Contig found in output but not in chromosomes : ', supercontig.names)
            return False
        
    #Then check if all true links still exist (i.e. links were not accidentally deleted)

    expectedContacts = sparse.dok_matrix((len(names), len(names))) #this is the matrix of all the links that are in the solution and therefore are expected in the output
    for c in chromosomes :
        for contig in range(len(c)-1) :
            expectedContacts[names[c[contig]],names[c[contig+1]]] = 1
            expectedContacts[names[c[contig+1]],names[c[contig]]] = 1
 
        #First, look if the expected links are found between supercontigs
    for segment in listOfSuperContigs :
        for endOfSegment in range(2):
            for l, neighbor in enumerate(segment.links[endOfSegment]) :
                expectedContacts[ names[segment.names[-endOfSegment]] , names[neighbor.names[-segment.otherEndOfLinks[endOfSegment][l]]]] = 0

        #Then within supercontigs
    for segment in listOfSuperContigs :
        for j in range(len(segment.listOfContigs)-1) :
            expectedContacts[names[segment.names[j]] , names[segment.names[j+1]]] = 0
            expectedContacts[names[segment.names[j+1]] , names[segment.names[j]]] = 0
            
    for i in range(len(names)):
        for j in range(len(names)):
            if expectedContacts[i,j] == 1 :
                print('In the actual chromosomes, there exist a link between ', i, ' and ', j, ' that is not found in the output')
                print(names)
                for k in listOfSuperContigs :
                    k.print_complete()
                return False
    
    return True


def stats_on_solve_ambiguities(n = 100, lengthOfChromosomes = 10, steps = 10) :
    
    record = []
    for i in range(n):
        chromosomes = buildFakeChromosomes(lengthOfChromosomes)
        lengthOfContig = 10000
        exportFakeToGFA(chromosomes, 'tests/stats/test' + str(i)+'.gfa', lengthOfContig)
        bf.export_to_csv(chromosomes, 'tests/stats/test' + str(i)+'.chro')
        
        listOfSegments, names = load_gfa('tests/stats/test' + str(i)+'.gfa')

        interactionMatrix = constructFakeInteractionMatrix(chromosomes, names, lengthOfContig)

        listOfSegments = solve_ambiguities(listOfSegments, interactionMatrix, 0.2, 0.45 ,5) #rejectedThreshold<AcceptedThreshold

        #record.append(check_result(chromosomes, listOfSuperContigs, names, links))
        
        bf.export_to_GFA(listOfSegments, exportFile = 'tests/stats/test' + str(i)+'F.gfa')
        #draw_distance_HiCcontacts_correlation(listOfSuperContigs, links, [10000 for i in names], interactionMatrix)
        
    fileRecord = open('tests/stats/record.txt', 'w')
    for i in record :
        fileRecord.write(str(i)+'\n')
    print(int(record.count(False)*100/n), '% of incorrectly changed GFA')

# chromosomes = ['A0-A1-A2-A3-A4-A5-A6-A7-A8-A9'.split('-'), 'A0-A1-A2-A3*-A4-A5-A6-A7-A8-A9'.split('-'),\
#                 'B0*-B1-B1-B2-B3-B4*-B5-B6-B7-B8-B9'.split('-'), 'B0*-B1-B2*-B3-B4-B5-B6-B7-B8-B9'.split('-')]

#chromosomes = bf.import_from_csv('tests/fake.chro')
# chromosomes = buildFakeChromosomes(10)
# #bf.export_to_csv(chromosomes, 'tests/fake.chro')

# lengthOfContig = 10000
# exportFakeToGFA(chromosomes, 'tests/fake.gfa', lengthOfContig)
# bf.export_to_csv(chromosomes, 'tests/fake.chro')
# listOfSegments, names = load_gfa('tests/fake.gfa')

# interactionMatrix = constructFakeInteractionMatrix(chromosomes, names, lengthOfContig)
# listOfSegments = solve_ambiguities(listOfSegments, interactionMatrix , 0.2, 0.45 ,5) #rejectedThreshold<AcceptedThreshold

# #flatten_loop(links, listOfSuperContigs, 1, 2)
# export_to_GFA(listOfSegments, gfaFile = 'tests/fake.gfa', exportFile = 'tests/fakeF.gfa', useExistingOffsetsFile = False)

# links, listOfSuperContigs, cn = simulated_annealing(originalLinks, names, interactionMatrix, [lengthOfContig for i in names], lambda x:1, 0.2, 0.45 ,5)
# export_to_GFA(links, listOfSuperContigs, cn, originalLinks, names = names, exportFile = 'tests/fakeA.gfa')

# print('Before the beginning of the process, the gfa energy is : ', \
#       score_output([[i] for i in range(len(names))], originalLinks, [lengthOfContig for i in names], interactionMatrix, infinite_distance = 500000))

#passing the dist_law is very inefficient, much too much redundant integration of this fucntion (gain of time possible)

#print('And the output is : ', check_result(chromosomes, listOfSegments, names))#, ', of energy ', score_output(listOfSuperContigs, links, [lengthOfContig for i in names], interactionMatrix, infinite_distance = 500000))

t = time.time()
stats_on_solve_ambiguities(n=200)

print('Finished in ', time.time()-t, ' seconds')
    