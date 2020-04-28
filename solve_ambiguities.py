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
from shutil import copyfile

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
def solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs):
            
    #we add all the new supercontigs
    for i in links[endOfSuperContig]:
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
    
       
    deletedContigs = [int(endOfSuperContig/2)]
    otherEnd = endOfSuperContig + 1- 2*(endOfSuperContig%2)
    for i in links[otherEnd]:
        links[i].remove(otherEnd)
    
    for merged in links[endOfSuperContig]:       
        if len(links[merged]) == 1 : #then the original copy is fully integrated in the supercontig
            print(endOfSuperContig, merged, links[merged+1-2*(merged%2)], links[merged])
            deletedContigs += [int(merged/2)]
            otherEnd = merged + 1- 2*(merged%2)
            for i in links[otherEnd]:
                links[i].remove(otherEnd)
        else : #then the original contig still exists by itself, we just delete the link
            links[merged].remove(endOfSuperContig)
    
    #then we replace the merged supercontigs and all their links by empty lists (we do not delete them to keep the indexes right)

    for i in deletedContigs :
        listOfSuperContigs[i] = [-1]
        links[2*i] = [-1]
        links[2*i+1] = [-1]
    
    return links, listOfSuperContigs

def merge_simply_two_adjacent_contig(endOfSuperContig, links, listOfSuperContigs):
    
    #we add the new supercontigs
    deleted = links[endOfSuperContig][0]
    listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig/2)]+listOfSuperContigs[int(deleted/2)]]
    
    otherEnd = endOfSuperContig +1 - 2*(endOfSuperContig%2)
    links += [links[otherEnd]]

    for j in links[otherEnd] :
        links[j] += [len(links)-1]
        
    otherEnd = deleted + 1 - 2*(deleted%2)
    links += [links[otherEnd]]
    for j in links[otherEnd] :
        links[j] += [len(links)-1]
            
    #now we delete the merged supercontig    
    deletedContigs = [int(endOfSuperContig/2), int(deleted/2)]
    otherEnd = endOfSuperContig + 1- 2*(endOfSuperContig%2)
    for i in links[otherEnd]:
        links[i].remove(otherEnd)
         
    otherEnd = deleted + 1- 2*(deleted%2)
    for i in links[otherEnd]:
        links[i].remove(otherEnd)
    
    #then we replace the merged supercontigs and all their links by empty lists (we do not delete them to keep the indexes right)
    for i in deletedContigs :
        listOfSuperContigs[i] = [-1]
        links[2*i] = [-1]
        links[2*i+1] = [-1]
        
    return links, listOfSuperContigs
    

def get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix) :
    
    for endOfContig in range(len(links)):
        if len(links[endOfContig]) > 1 :
            absoluteLinksStrength, linksStrength = intensity_of_interactions(listOfSuperContigs[int(endOfContig/2)], \
                                                                             [listOfSuperContigs[int(x/2)] for x in links[endOfContig]], interactionMatrix)
            maxStrength = np.max(linksStrength)
            
            for i in range(len(links[endOfContig])-1,-1,-1) :
                if linksStrength[i] < maxStrength/3 : #we consider then that the link does not exist
                    links[links[endOfContig][i]].remove(endOfContig)
                    del links[endOfContig][i]
    return links

def clean_listOfSuperContigs(links, listOfSuperContigs):
     
    if links[0] == [-1]:
        list_of_number_of_bad_lists= [1]
    else :
        list_of_number_of_bad_lists= [0]
    
    for i in links[1:] :
        if i == [-1] :
            list_of_number_of_bad_lists += [list_of_number_of_bad_lists[-1]+1]
        else :
            list_of_number_of_bad_lists += [list_of_number_of_bad_lists[-1]]
            
    #now we get rid of all bad lines
    links = [i for i in links if i!=[-1]]
    listOfSuperContigs = [i for i in listOfSuperContigs if i!=[-1]]
    
    for i in range(len(links)) :
        for j in range(len(links[i])) :
            links[i][j] = links[i][j] - list_of_number_of_bad_lists[links[i][j]]
            
    return links, listOfSuperContigs
    
def merge_contigs (links, listOfSuperContigs):
    #each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends)
    locked = [False for i in listOfSuperContigs]
    
    for endOfSuperContig in range(len(locked)*2):
        #print('end of supercontig that is looked at : ', endOfSuperContig)
        if len(links[endOfSuperContig]) > 1 :
            startMerging = not locked[int(endOfSuperContig/2)]
            for i in links[endOfSuperContig] :
                if i < len(locked)*2 :
                    startMerging = startMerging and (not(locked[int(i/2)]))
                else :
                    startMerging = False
            
            if startMerging : #if nothing is locked for now
                links, listOfSuperContigs = solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs)
                links = [[links[i][j] for j in range(len(links[i]))] for i in range(len(links))] #we might want to see where that comes from if we want speed
                locked[int(endOfSuperContig/2)] = True
                for i in links[endOfSuperContig] :
                    locked[int(i/2)] = True
            
    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    
    for endOfSuperContig in range(len(locked)*2):
        if len(links[endOfSuperContig]) == 1 and len(links[links[endOfSuperContig][0]])==1 : #then we just merge
            links, listOfSuperContigs = merge_simply_two_adjacent_contig(endOfSuperContig, links, listOfSuperContigs)
        
    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    return links, listOfSuperContigs
            
def export_to_GFA(links):
    
    copyfile('results/sequencesWithoutLinks.gfa', 'results/newAssembly.gfa')
    f = open('results/newAssembly.gfa', 'a')
    
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
    steps = 1
    for i in range(steps) :
        if i%100 == 0 :
            print (str(i/steps*100)+'% done')
        
        links = get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix)
        print('Bad links deleted')
        links, listOfSuperContigs = merge_contigs(links, listOfSuperContigs)
        
    return links, listOfSuperContigs 

links = bf.import_links('listsPython/links.csv')
print(links[1465])
#infContigs = bf.read_info_contig('data/results/info_contigs.txt')
interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')
for i in range(len(interactionMatrix)):
    interactionMatrix[i][i] = 0
print('Loaded')

#links, listOfSuperContigs = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix)
#print(listOfSuperContigs)
#print('Now the end of links ')
#print(links[2600:])
#print(intensity_of_interactions([874*2, 874*2+1], [[584*2,584*2+1,1120*2,1120*2+1], [584*2,584*2+1,78*2,78*2+1]], interactionMatrix))
#print(intensity_of_interactions([584*2,584*2+1,78*2,78*2+1], [[802*2,802*2+1,874*2, 874*2+1], [802*2,802*2+1,743*2,743*2+1]], interactionMatrix))

#print(solve_ambiguities([[],[4],[],[4],[1,3],[6,8],[5],[],[5],[]], [0,1,2,3,4], [[1,1,1,1,0],[1,1,1,0,1],\
#                                                                                 [1,1,1,1,1],[1,0,1,1,1],\
#                                                                                     [0,1,1,1,1]]))

links, listOfContigs = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix)
print(listOfContigs)    
                         
print('Finished')

    