#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from transform_gfa import load_gfa
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

def HiC_vs_GFA(hiccontacts, links, fragment_list):

    confirmationOfLinks = [[0 for i in j] for j in links]  # a list of list of 0 of the same dimensions as links

    for contact in hiccontacts:
        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]

        for j in range(len(links[contig1 * 2])):
            if (
                links[contig1 * 2][j] == contig2 * 2
                or links[contig1 * 2][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2][j] += contact[2]
                for i in range(len(links[links[contig1 * 2][j]])):
                    if links[links[contig1 * 2][j]][i] == contig1 * 2:
                        confirmationOfLinks[links[contig1 * 2][j]][i] += contact[2]

        for j in range(len(links[contig1 * 2 + 1])):
            if (
                links[contig1 * 2 + 1][j] == contig2 * 2
                or links[contig1 * 2 + 1][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                for i in range(len(links[links[contig1 * 2 + 1][j]])):
                    if links[links[contig1 * 2 + 1][j]][i] == contig1 * 2 + 1:
                        confirmationOfLinks[links[contig1 * 2 + 1][j]][i] += contact[2]

    return confirmationOfLinks

# same as above but taking into account contigs that are two connexions away
def HiC_vs_GFAtwo(hiccontacts, links, fragment_list, coverage):  
    
    confirmationOfLinks = [
        [0 for i in j] for j in links
    ]  # a list of list of 0 of the same dimensions as links
    weightedconfirmationOfLinks = [[0 for i in j] for j in links]

    for contact in hiccontacts:

        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]

        for j, neighbor in enumerate(links[contig1 * 2]):
            # direct neighbor
            if (
                links[contig1 * 2][j] == contig2 * 2
                or links[contig1 * 2][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2][j] += contact[2]
                weightedconfirmationOfLinks[contig1 * 2][j] += (
                    contact[2] / coverage[contig1] / coverage[contig2]
                )

                for i in range(len(links[links[contig1 * 2][j]])):
                    if links[links[contig1 * 2][j]][i] == contig1 * 2:
                        confirmationOfLinks[links[contig1 * 2][j]][i] += contact[2]
                        weightedconfirmationOfLinks[links[contig1 * 2][j]][i] += (
                            contact[2] / coverage[contig1] / coverage[contig2]
                        )
            # two connexions away
            for c in range(
                len(links[neighbor + 1 - 2 * neighbor % 2])
            ):  # we take the other end of the neighbor contig
                if (
                    links[neighbor + 1 - 2 * neighbor % 2][c] == contig2 * 2
                    or links[neighbor + 1 - 2 * neighbor % 2][c] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    confirmationOfLinks[neighbor + 1 - 2 * neighbor % 2][c] += contact[
                        2
                    ]
                    weightedconfirmationOfLinks[neighbor + 1 - 2 * neighbor % 2][c] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(
                        len(links[links[neighbor + 1 - 2 * neighbor % 2][c]])
                    ):
                        if (
                            links[links[neighbor + 1 - 2 * neighbor % 2][c]][i]
                            == neighbor + 1 - 2 * neighbor % 2
                        ):
                            confirmationOfLinks[
                                links[neighbor + 1 - 2 * neighbor % 2][c]
                            ][i] += contact[2]
                            weightedconfirmationOfLinks[
                                links[neighbor + 1 - 2 * neighbor % 2][c]
                            ][i] += (contact[2] / coverage[contig1] / coverage[contig2])
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1 * 2:
                            confirmationOfLinks[neighbor][i] += contact[2]
                            weightedconfirmationOfLinks[neighbor][i] += (
                                contact[2] / coverage[contig1] / coverage[contig2]
                            )

        for j, neighbor in enumerate(links[contig1 * 2 + 1]):
            for j in range(len(links[contig1 * 2 + 1])):
                if (
                    links[contig1 * 2 + 1][j] == contig2 * 2
                    or links[contig1 * 2 + 1][j] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2 + 1][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(len(links[links[contig1 * 2 + 1][j]])):
                        if links[links[contig1 * 2 + 1][j]][i] == contig1 * 2 + 1:
                            confirmationOfLinks[links[contig1 * 2 + 1][j]][
                                i
                            ] += contact[2]
                            weightedconfirmationOfLinks[links[contig1 * 2 + 1][j]][
                                i
                            ] += (contact[2] / coverage[contig1] / coverage[contig2])

            # two connexions away
            otherEndOfNeighbor = neighbor + 1 - 2 * (neighbor % 2)

            for c in range(
                len(links[otherEndOfNeighbor])
            ):  # we take the other end of the neighbor contig
                if (
                    links[otherEndOfNeighbor][c] == contig2 * 2
                    or links[otherEndOfNeighbor][c] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2 + 1][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    confirmationOfLinks[otherEndOfNeighbor][c] += contact[2]
                    weightedconfirmationOfLinks[otherEndOfNeighbor][c] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(len(links[links[otherEndOfNeighbor][c]])):

                        if links[links[otherEndOfNeighbor][c]][i] == otherEndOfNeighbor:
                            confirmationOfLinks[links[otherEndOfNeighbor][c]][
                                i
                            ] += contact[2]
                            weightedconfirmationOfLinks[links[otherEndOfNeighbor][c]][
                                i
                            ] += (contact[2] / coverage[contig1] / coverage[contig2])
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1 * 2:
                            confirmationOfLinks[neighbor][i] += contact[2]
                            weightedconfirmationOfLinks[neighbor][i] += (
                                contact[2] / coverage[contig1] / coverage[contig2]
                            )
                    # if verif != 2 :
                    #   print ('il y a un probleme')

    return confirmationOfLinks, weightedconfirmationOfLinks

# a function to delete small contigs made of repeated sequences that have no HiC contacts but tons of links
def crush_small_contigs(segments, interactionMatrix) :
    
    #first list all segments to be deleted
    for segment in segments :
        if segment.length < 5000 and segment in segment.links[0] : #this characterize a small repeated sequence
            if segment.HiCcoverage == 0 : # this segment is blind, can't do anything but crush it
                segment.ID = -segment.ID #segments marked with negative IDs should be deleted
                print(segment.names, ' is a small contig')
                
    print('Listed all segments to delete')
            
    #then delete all the links going toward unwanted segments (we don't want to reroute useless links)
    for segment in segments :
        for endOfSegment in range(2) :
            for n in range(len(segment.links[endOfSegment])-1, -1, -1) :
                if segment.links[endOfSegment][n].ID < 0 :
                    del segment.links[endOfSegment][n]
                    del segment.otherEndOfLinks[endOfSegment][n]
                    del segment.CIGARs[endOfSegment][n]                                        
    print('Deleted all links towards bad segments')
    
    #reroute all links going in and out of crushed segments
    # for segment in segments :
    #     if segment.ID < 0 :
    
    #         print('Crushing segment ', segment.names)
            
    #         #then make this segment disappear, by linking segments to the left of this contig to segment to the right :
    #         for l, leftneighbor in enumerate(segment.links[0]) :
    #             if leftneighbor.ID != segment.ID :
    #                 #print('Ln before : ', leftneighbor.ID, [se.ID for se in leftneighbor.links[segment.otherEndOfLinks[0][l]]])
    #                 leftneighbor.add_a_bunch_of_end_of_links(segment.otherEndOfLinks[0][l], segment.links[1], segment.otherEndOfLinks[1], ['0M' for i in segment.links[1]])
                    
    #         for r, rightneighbor in enumerate(segment.links[1]):
    #             if rightneighbor.ID != segment.ID :
    #                 #print('Rn before : ', [se.ID for se in rightneighbor.links[segment.otherEndOfLinks[1][l]]])
    #                 rightneighbor.add_a_bunch_of_end_of_links(segment.otherEndOfLinks[1][r], segment.links[0], segment.otherEndOfLinks[0], ['0M' for i in segment.links[0]])
    
    # print('Rerouted all links around segments to delete')
    #delete all segments that should be
    for se in range(len(segments)-1, -1, -1) :
        if segments[se].ID < 0 :
            del segments[se]
            
    print('There are ', len(segments), ' segments left')
                     
                

