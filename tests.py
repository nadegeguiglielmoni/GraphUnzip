#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:05:13 2020

@author: zaltabar

In this file, test functions to test our algorithm
"""

import matplotlib.pyplot as plt
import random
import numpy as np

import basic_functions as bf
from transform_gfa import load_gfa
from basic_functions import export_to_GFA
from solve_ambiguities import solve_ambiguities
from solve_ambiguities import intensity_of_interactions
import scipy.integrate as integrate

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

def exportFakeToGFA(chromosomes, file) :
    
    f = open(file, 'w')
    
    alreadyExported = []
    for c in chromosomes :
        for contig in c :
            if contig not in alreadyExported :
                f.write('S\t'+contig+'\t*\n')
                alreadyExported += [contig]
      
    linksAlreadyExported = []
    for c in chromosomes :
        for contig in range(len(c)-1) :
            if [c[contig],c[contig+1]] not in linksAlreadyExported :
                linksAlreadyExported += [[c[contig],c[contig+1]]]
                f.write('L\t'+c[contig]+'\t+\t'+c[contig+1]+'\t+\t*\n')
                
    return alreadyExported
      
def dist_law(distance) :
    if distance < 10000 :
        distance = 5000
    return 10000/distance

def constructFakeInteractionMatrix(chromosomes, names, lengthOfContigs = 10000):
    
    interactionMatrix = [[0 for i in names] for j in names]
    for c in chromosomes :
        for c1 in range(len(c)):
            for c2 in range(c1, len(c)) :
                con1 = c[c1]
                con2 = c[c2]
                interactionMatrix[names.index(con1)][names.index(con2)] += integrate.quad(dist_law, (c2-c1-1)*lengthOfContigs, (c2-c1)*lengthOfContigs)[0]
                interactionMatrix[names.index(con2)][names.index(con1)] += integrate.quad(dist_law, (c2-c1-1)*lengthOfContigs, (c2-c1)*lengthOfContigs)[0]
    
    for i in range(len(interactionMatrix)):
        interactionMatrix[i][i] = 0
    return interactionMatrix

#function that checks if ls1 is a sublist of ls2 : useful for checking the quality of the output
def sublist(ls1, ls2):

    s1 = ','.join(ls1)
    s2 = ','.join(ls2)

    return s1 in s2

#function taking as arguments the solution of the problem and the output of the algorithm to see if the output is wrong
#for now, it just checks if all supercontig in the output actually exist, but some links may be missing
def check_result(chromosomes, listOfSuperContigs, names) :
    for supercontig in listOfSuperContigs :
        
        supercontigname = [names[i] for i in supercontig]
        found = False
        for c in chromosomes :
            if sublist(supercontigname,c) :
                found = True
        if not found :
            for c in chromosomes :
                print(c)
            print('Contig found in output but not in chromosomes : ', supercontigname)
            return False
    return True
        
# chromosomes = buildFakeChromosomes(10)

chromosomes = ['A0*-A1-A2-A3-A4-A5-A6*-A7-A8-A9'.split('-'), 'A0-A1-A2-A3-A4-A5*-A6-A7-A8-A9*'.split('-'),\
                'B0-B1-B2-B3-B4-B5-B6-B7-B8-B9'.split('-'), 'B0-B1-B2-B0-B3-B4-B5-B6-B7-B8-B9*'.split('-')]

exportFakeToGFA(chromosomes, 'tests/fake.gfa')
bf.export_to_csv(chromosomes, 'tests/fake.chro')

#chromosomes = bf.import_from_csv('tests/historyOfCases/failure6.chro')

links, names = load_gfa('tests/fake.gfa')

print(chromosomes)


lengthOfContig = 10000
interactionMatrix = constructFakeInteractionMatrix(chromosomes, names, lengthOfContig)
print_chromosomes(chromosomes)
print(names)

links, listOfSuperContigs, cn = solve_ambiguities(links, names, interactionMatrix, [lengthOfContig for i in names], dist_law, 0.2, 0.45 ,10) #rejectedThreshold<AcceptedThreshold
#passing the dist_law is very inefficient, much too much redundant integration of this fucntion (gain of time possible)

print('And the output is : ', check_result(chromosomes, listOfSuperContigs, names))


#export_to_GFA(links, listOfSuperContigs, cn, names, exportFile = 'tests/fake2.gfa')
#export_to_GFA(links, [[i] for i in range(len(names))], [1 for i in names], names, exportFile = 'tests/reexport.gfa')
print('Finished')
    