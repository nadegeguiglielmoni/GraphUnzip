#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

@author: zaltabar

In this file we will use a lot of "supercontigs", meaning an assembly of several contigs. The syntax is a list
of end of contigs, e.g. [40,41, 104,105,523,522] (and not [20,52,261])
"""

import matplotlib.pyplot as plt
import numpy as np
import basic_functions as bf

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
                if c%2 == 0 and contig%2 == 0  :#to count each interaction only once
                    if c not in commonContigs :
                        absoluteScore[sg] += interactionMatrix[int(c/2)][int(contig/2)]
                        relativeScore[sg] += interactionMatrix[int(c/2)][int(contig/2)]
                    else :
                        absoluteScore[sg] += interactionMatrix[int(c/2)][int(contig/2)]
    
    return absoluteScore, relativeScore
    
#we're going to look specifically at one contig and its immediate surroundings
def solve_ambiguity_around_this_contig(contig, links, interactionMatrix):
         
    linksStrengthEvenEnd = intensity_of_interactions([contig*2, contig*2+1], [[x, x +1-2*(x%2)] for x in links[contig*2]], interactionMatrix)[1]
    maxStrength = np.max(linksStrengthEvenEnd)
    
    newSuperContigs = []
    for i in range(len(linksStrengthEvenEnd)) :
        if linksStrengthEvenEnd[i] > maxStrength/4 :
            newSuperContigs += [[links[contig*2][i], links[contig*2][i] +1-2*(links[contig*2][i]%2), contig*2, contig*2+1]]
    print (newSuperContigs)

links = bf.import_links('listsPython/links.csv')
#infContigs = bf.read_info_contig('data/results/info_contigs.txt')
interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')

print('coucou')
solve_ambiguity_around_this_contig(802, links, interactionMatrix)
#print(intensity_of_interactions([874*2, 874*2+1], [[584*2,584*2+1,1120*2,1120*2+1], [584*2,584*2+1,78*2,78*2+1]], interactionMatrix))
#print(intensity_of_interactions([584*2,584*2+1,78*2,78*2+1], [[802*2,802*2+1,874*2, 874*2+1], [802*2,802*2+1,743*2,743*2+1]], interactionMatrix))

print('Finished')

    