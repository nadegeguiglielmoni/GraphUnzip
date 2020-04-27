#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

@author: zaltabar
"""

import matplotlib.pyplot as plt
import numpy as np
import basic_functions as bf
import random

def how_far_away_are_those_contigs(contig1, contig2, links, infContigs):
    connectedToContig1 = [contig1*2, contig1*2+1]
    
    distanceToContig1 = [0,0]
    
    pathToContig1 = [[],[]]
    
    pastLengthOfconnexion1 = 1
    
    while pastLengthOfconnexion1 < len(connectedToContig1) :
        
        pastLengthOfconnexion1 = len(connectedToContig1)
        newConnectedToContig1 = connectedToContig1[:]
        
        #we spread one step from contig1
        for cc, connectedContig in enumerate(connectedToContig1) : #that can be optimized if we need speed, since we go through the first contigs many times
            
            if cc < 2 or cc%2 == 1 : #that's to move only one way in the contig graph
                for newContig in links[connectedContig] :
                            
                    newdist = distanceToContig1[cc]+infContigs[int(connectedContig/2)][1]
                    if newContig in newConnectedToContig1 :
                        
                        oldpos = newConnectedToContig1.index(newContig)
                        
                        if newdist < distanceToContig1[oldpos] :
                            distanceToContig1[oldpos] = newdist
                            pathToContig1[oldpos] = pathToContig1[cc]+[cc]
                            
                            distanceToContig1[oldpos+1-2*(oldpos%2)] = newdist
                            pathToContig1[oldpos+1-2*(oldpos%2)] = pathToContig1[cc]+[cc]
                    else :
                        #we add both end of the contig to the connected contig
                        newConnectedToContig1 += [newContig]
                        distanceToContig1 += [newdist]
                        pathToContig1 += [pathToContig1[cc]+[cc]]
                        
                        newConnectedToContig1 += [newContig + 1 - 2*(newContig%2)]
                        distanceToContig1 += [newdist]
                        pathToContig1 += [pathToContig1[cc]+[connectedContig]]
        connectedToContig1 = newConnectedToContig1[:]
        
    #now we see where contig2 stands with respect to contig1
    if contig2*2 in connectedToContig1 :
        print(connectedToContig1)
        indexContig2 = connectedToContig1.index(contig2*2)
        return distanceToContig1[indexContig2]-infContigs[contig1][1], pathToContig1[indexContig2][1:] #because the length of contig 1 is always taken into account in the path length while it shouldn't
    else :
        print(connectedToContig1)
        return -1 #meaning contig1 and contig2 are not connected in this graph

#this function measures the intensity of interactions between one supercontig and several candidate, including without taking account of the common parts of the supercontigs
#the contig in this function are not numbered by their end, i.e. give it [1234] and not [2468,2469]
def intensity_of_interactions(supercontig, listOfSuperContigs, interactionMatrix) :
    commonContigs = []
    for contig in listOfSuperContigs[0] :
            common = True
            for sg in listOfSuperContigs[1:] :
                common = common and (contig in sg)
            if common :
                commonContigs += [contig]
    #now we have a list of common contigs
    
    #we return for each supercontig its absolute score and its relative score (wihtout the common parts)
    absoluteScore = [0 for i in range(len(listOfSuperContigs))]
    relativeScore = [0 for i in range(len(listOfSuperContigs))]
    
    for sg in range(len(listOfSuperContigs)) :
        for c in listOfSuperContigs[sg] :
            for contig in supercontig :
                if c not in commonContigs :
                    absoluteScore[sg] += interactionMatrix[c][contig]
                    relativeScore[sg] += interactionMatrix[c][contig]
                else :
                    absoluteScore[sg] += interactionMatrix[c][contig]
    
    return absoluteScore, relativeScore
    
#we're going to look specifically at one contig and its immediate surroundings
def solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs, interactionMatrix):
         
    absoluteLinksStrength, linksStrength = intensity_of_interactions([int(endOfSuperContig/2)], [[int(x/2)] for x in links[endOfSuperContig]], interactionMatrix)
    maxStrength = np.max(linksStrength)
    
    validatedLinks = []
    for i in range(len(linksStrength)) :
        if linksStrength[i] >= maxStrength/3 : #we consider then that the link is real
            validatedLinks += [links[endOfSuperContig][i]]
            print('Validated : ' , links[endOfSuperContig][i])
        else : #we consider that the link does not exist, and remove one end, we'll remove the other one after
            links[links[endOfSuperContig][i]].remove(endOfSuperContig)
    #removing the other end of rejected links
    for i in range(len(links[endOfSuperContig])-1,-1,-1):
            if links[endOfSuperContig][i] not in validatedLinks :
                del links[endOfSuperContig][i]
            
    #we add all the new supercontigs
    for i in validatedLinks:
        listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig/2)]+listOfSuperContigs[int(i/2)]]
        
        otherEnd = endOfSuperContig +1 - 2*(endOfSuperContig%2)
        links += [links[otherEnd]]
        for j in links[otherEnd] :
            links[j] += [len(links)-1]
            
        otherEnd = i + 1 - 2*(i%2)
        links += [links[otherEnd]]
        for j in links[otherEnd] :
            links[j] += [len(links)-1]
            
    #now we delete the merged supercontigs
    #we start by deleting the links that linked the merged supercontigs to the outside
    otherEnd = endOfSuperContig+1-2*(endOfSuperContig%2)
    for oldLink in range(len(links[otherEnd])):
            links[links[otherEnd][oldLink]].remove(otherEnd)
    for merged in validatedLinks:
        otherEnd = merged+1-2*(merged%2)
        for oldLinks in range(len(links[otherEnd])):
            links[links[otherEnd][oldLinks]].remove(otherEnd)
    
    #then we replace the merged supercontigs and all their link by empty lists (we do not delete them to keep the indexes right)
    deletedContigs = [int(endOfSuperContig/2)]+[int(i/2) for i in validatedLinks]

    for i in range(len(links)) :
        if int(i/2) in deletedContigs:
            links[i] = []
    for i in range(len(listOfSuperContigs)):
        if i in deletedContigs :
            listOfSuperContigs[i] = []
    
    return links, listOfSuperContigs

def export_to_GFA(links, listOfContigs, fastaFile):
    
    f = open('results/newAssembly.gfa', 'w')
    f.write('H\tVN:Z:1.0\n')
    
    for i in range(len(listOfContigs)) :
        if i%200 == 0 :
            print(i)
        f.write('S\t'+str(i*2)+'\t'+bf.get_contig(fastaFile, listOfContigs[i])+'\tR:i:'+str(i*10000)+'\n')
    for i in range(len(links)) :
        for j in range(len(links[i])):
            if i < links[i][j] :
                if i%2 == 0 and links[i][j]%2 == 0:
                    f.write('L\t'+str(i)+'\t-\t'+str(links[i][j])+'\t+\t*\n')
                elif i%2 == 1 and links[i][j]%2 == 0:
                    f.write('L\t'+str(i-1)+'\t+\t'+str(links[i][j])+'\t+\t*\n')
                elif i%2 == 0 and links[i][j]%2 == 1:
                    f.write('L\t'+str(i)+'\t-\t'+str(links[i][j]-1)+'\t-\t*\n')
                elif i%2 == 1 and links[i][j]%2 == 1:
                    f.write('L\t'+str(i-1)+'\t+\t'+str(links[i][j]-1)+'\t-\t*\n')

def solve_ambiguities(links, listOfContigs, interactionMatrix): #look at ambilguities one after the other 
    listOfSuperContigs = [[x] for x in listOfContigs]
    steps = 5
    for i in range(steps) :
        if i%100 == 0 :
            print (str(i/steps*100)+'% done')
        inspectedContig = random.randint(0, len(links)-1)
        if len(links[inspectedContig]) > 1 :
            print('Let us try and solve contig number '+str(inspectedContig))
            links, listOfSuperContigs = solve_ambiguity_around_this_end_of_contig(inspectedContig, links, listOfSuperContigs, interactionMatrix)
    
    return links, listOfSuperContigs 

links = bf.import_links('listsPython/links.csv')
#infContigs = bf.read_info_contig('data/results/info_contigs.txt')
interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')
print('Loaded')

links, listOfSuperContigs = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix)
print(listOfSuperContigs)
print('Now the end of links ')
print(links[2600:])
#print(intensity_of_interactions([874*2, 874*2+1], [[584*2,584*2+1,1120*2,1120*2+1], [584*2,584*2+1,78*2,78*2+1]], interactionMatrix))
#print(intensity_of_interactions([584*2,584*2+1,78*2,78*2+1], [[802*2,802*2+1,874*2, 874*2+1], [802*2,802*2+1,743*2,743*2+1]], interactionMatrix))
#export_to_GFA(links, listOfContigs, 'data/Assembly.fasta')

print('Finished')

    