#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:22:21 2020

@author: zaltabar

This file is for function that do not have anymore utility in the master code, but that took time to program
"""

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

def detect_fishy_links(links, confirmationOfLinks, coverage):

    # we're going to detect, when there is an ambiguity, if one path looks unlikely
    badlinks = [] * len(links)
    for i in coverage:  # to ensure we don't get absurdly high multiplicative factors
        if i < 0.01:
            i = 0.01

    for endOfContig in range(len(links)):

        if len(links[endOfContig]) > 1:

            weightedConfirmation = [-1] * len(links[endOfContig])
            for i in range(len(links[endOfContig])):

                if (
                    coverage[int(endOfContig / 2)] > 0.01
                    and coverage[int(links[endOfContig][i] / 2)] > 0.01
                ):
                    weightedConfirmation[i] = (
                        confirmationOfLinks[endOfContig][i]
                        / coverage[int(links[endOfContig][i] / 2)]
                    )

            maximum = np.max(weightedConfirmation)
            reliable = -1 not in weightedConfirmation

            if reliable:
                for i in range(len(weightedConfirmation)):
                    if weightedConfirmation[i] < 0.1 * maximum:

                        badlinks += [[endOfContig, i]]
                        print(
                            "There is a suspect link here : ",
                            endOfContig,
                            links[endOfContig],
                            confirmationOfLinks[endOfContig],
                            weightedConfirmation,
                        )
    return badlinks