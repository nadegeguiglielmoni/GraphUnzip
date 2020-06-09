#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import numpy as np
import basic_functions as bf
import scipy.integrate as integrate
from bisect import bisect_left #to look through sorted lists

from transform_gfa import check_links
import segment as s
from segment import Segment

# this function measures the intensity of interactions between one supercontig and several candidate, including without taking account of the common parts of the supercontigs
# It also weighs the interaction with the length of a supercontig, so that a very long candidate spercontig is not seen as having a lot of connexion just because it is long
def intensity_of_interactions(
    segment,
    candidatesSegments,
    listOfSegments,
    interactionMatrix,
    dist_law,
    copiesnumber,
    supercontigsaretouching=False,
):
    for candidate in candidatesSegments :
        if candidate == segment : #small loop, don't solve that !
            return [-1], [-1]
    
    commonContigs = []
    potentialCommonContigs = candidatesSegments[0].listOfContigs.copy()
    if supercontigsaretouching:
        if segment in candidatesSegments[0].links[0] :
            for i in candidatesSegments[0].links[1] :
                potentialCommonContigs += i.listOfContigs
        elif segment in candidatesSegments[0].links[1] :
            for i in candidatesSegments[0].links[0] :
                potentialCommonContigs += i.listOfContigs
        else :
            supercontigsaretouching = False
            print('WARNING in intensity_of_interactions : called like the contigs were touching but they were not')
            candidatesSegments[0].print_complete()
            segment.print_complete()

    for contig in potentialCommonContigs:
        
        common = True

        for n, candidate in enumerate(candidatesSegments[1:]):

            allcontigs = candidate.listOfContigs.copy()
            if supercontigsaretouching:
                
                if segment in candidate.links[0] :
                    allcontigs += candidate.links[1].copy()
                elif segment in candidate.links[1] :
                    allcontigs += candidate.links[0].copy()
                else :
                    supercontigsaretouching = False
                    print('WARNING in intensity_of_interactions : called like the contigs were touching but they were not')
                    candidate.print_complete()
                    segment.print_complete()
            common = common and (contig in allcontigs)

        if common:
            commonContigs += [contig]

    # now we have a list of common contigs
    for candidate in candidatesSegments :
        if all(elem in commonContigs for elem in candidate.listOfContigs):  # if all elements of sg are in commoncontigs, the algorithm cannot make a choice for now
            return [-1], [-1]

    #bestsignature decides what contig is most characterisic of the segment
    bestSignature = np.min([copiesnumber[x] for x in segment.listOfContigs])
    # we return for each supercontig its absolute score and its relative score (wihtout the common parts)
    absoluteScores = []
    relativeScores = []

    for c in candidatesSegments:

        absoluteScore, relativeScore, partial_area = c.interaction_with_contigs(segment, interactionMatrix, dist_law, copiesnumber, commonContigs, bestSignature)

        absoluteScores.append(absoluteScore)
        relativeScores.append(relativeScore)
        
        if supercontigsaretouching:
            relativeScores[-1] /= partial_area

    #print("*Common contigs : ", commonContigs)
    return absoluteScores, relativeScores

def isPresent(l, x):
    i = bisect_left(l, x)
    if i != len(l) and l[i] == x:
        return True
    return False

# here we look specifically at one contig and its immediate surroundings (can return -1 if fails in short loop)
def duplicate_around_this_end_of_contig(segment, endOfSegment, listOfSuperContigs, copiesnumber): #endOfSegment should be 0 if it's the left end and1 if it's the right end
    
    if segment in segment.links[endOfSegment] : #if a segment loops on itself, another module would be needed
        return 0
    if any(segment.links[endOfSegment].count(i)>=2 for i in segment.links[endOfSegment]): #if segment.links[endOfSegment] has two copies of the same segment, it means one link going towards each end, that is not solvable
        return 0

    for i in segment.listOfContigs:
        copiesnumber[i] += len(segment.links[endOfSegment]) - 1

    # add all the new supercontigs
    for neighbor in segment.links[endOfSegment]:
        s.merge_two_segments(segment, endOfSegment, neighbor, listOfSuperContigs) #the merged segment is appended at the end of listOfSuperContigs

    # now delete the merged supercontigs
    # start by deleting the links that linked the merged supercontigs to the outside

    deletedContigs = [segment.hash()]
    otherEnd = 1 - endOfSegment

    for i, neighbor in enumerate(segment.links[otherEnd]):
    
        neighbor.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)
        

    for m, merged in enumerate(segment.links[endOfSegment]):
        
        if len(merged.links[segment.otherEndOfLinks[endOfSegment][m]]) == 1:  # then the original copy is fully integrated in the supercontig
            deletedContigs.append(merged.hash())
            
            otherEnd = 1-segment.otherEndOfLinks[endOfSegment][m]
            for i, neighbor in enumerate(merged.links[otherEnd]):
                try:
                    neighbor.remove_end_of_link(merged.otherEndOfLinks[otherEnd][i], merged, otherEnd)
                except ValueError:  # that means we're in a small loop which we can't solve
                    print("There is merging difficulty around the far end of "
                        + str(merged.names)+ " from "
                        + str(segment.names)+ " . Please check that there is indeed a loop there.")
                    return 0
                
        else:  # then the original contig still exists by itself, just delete the link going toward segment
            try:
                merged.remove_end_of_link(segment.otherEndOfLinks[endOfSegment][m], segment, endOfSegment)
            except ValueError:  # that means we're in a small loop which whe can't solve
                print("There is merging difficulty around the near end of "
                        + str(merged.names)+ " from "
                        + str(segment.names)+ " . Please check that there is indeed a loop there.")
                return 0

    #delete all segments that should be
    
    deletedContigs.sort()
    for i in range(len(listOfSuperContigs)-1,-1,-1) :
        h = listOfSuperContigs[i].hash()
        if isPresent(deletedContigs, h) : #in other words, if h is in deletedContigs (written like that because it has logarithmic efficiency) 
            del listOfSuperContigs[i]
            
    #lock all the segments that have been duplicated, so that they are not duplicated by both ends
    segment.lockNode(endOfSegment)

# similar to the function above, but simpler : put in one supercontig two smaller supercontig linked by a link unambinguous at both ends
def merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments):
  
    if len(segment.links[endOfSegment]) != 1 :
        print('ERROR : trying to merge simply two contigs that cannot be merged simply')
        return -1, -1
    
    neighbor = segment.links[endOfSegment][0]
    endOfSegmentNeighbor = segment.otherEndOfLinks[endOfSegment][0]

    if len(neighbor.links[endOfSegmentNeighbor]) != 1 :
        print('ERROR : trying to merge simply two contigs that cannot be merged simply')
        return -1,-1
        
    if neighbor == segment :  # then we do not merge a contig with itself
        return -1, -1

    # add the new segment
    s.merge_two_segments(segment, endOfSegment, neighbor, listOfSegments)

    # delete links going towards the two ex-segments
    otherEnd = 1 - endOfSegment
    otherEndNeighbor = 1 - endOfSegmentNeighbor
    
    for i, n in enumerate(segment.links[otherEnd]) :
        n.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment)
        
    for i, n in enumerate(neighbor.links[otherEndNeighbor]) :
        #print('Removing ', neighbor.names, ' from ', n.names, ' and adding the new contig',listOfSegments[-1].names, ' at end ', neighbor.otherEndOfLinks[otherEndNeighbor][i])
        n.remove_end_of_link(neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor)

    # delete the ex-segments
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)

    return listOfSegments

#a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):
        
    goOn = True
    while goOn :
        goOn = False
        for segment in listOfSegments:
            
            alreadyDidThisOne = False #if the segment is deleted when looking at its first end, you don't want it to look at its other end, since it does not exist anymore
            for endOfSegment in range(2):
                if not alreadyDidThisOne :
                    alreadyDidThisOne = True
                    if len(segment.links[endOfSegment]) == 1\
                        and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]]) == 1:  # then merge
                        if segment != segment.links[endOfSegment][0]:
                            goOn = True
                            listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)

    return listOfSegments

#get_rid_of_bad_links compare links using HiC contact informations when there is a choice and delete links that are not supported by HiC evidence
def get_rid_of_bad_links(listOfSegments,interactionMatrix,dist_law,copiesnumber,thresholdRejected,thresholdAccepted):

    # loop through all segments inspecting the robustness of all links.
    for segment in listOfSegments:
        
        for endOfSegment in range(2):
            
            if len(segment.links[endOfSegment]) >= 2 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
                
                # comparison pairwise of the links, those that should be deleted are deleted
                for n1 in range(len(segment.links[endOfSegment]) - 1):
                    n2 = n1 + 1
                    while n2 < len(segment.links[endOfSegment]):
                        absoluteLinksStrength, linksStrength = intensity_of_interactions(segment, [segment.links[endOfSegment][n1], segment.links[endOfSegment][n2]]\
                                                                                         , listOfSegments, interactionMatrix, dist_law, copiesnumber, True)
                        if absoluteLinksStrength == [-1]: #means that the configuration does not enable the algorithm to compare the two interactions
                            segment.freezeNode(endOfSegment)
                            
                        else:  
                            print('I have to decide, at ', segment.names, ' between ', segment.links[endOfSegment][n1].names, ' and ', segment.links[endOfSegment][n2].names, ' with these values : ', linksStrength)
                            if linksStrength[0] > linksStrength[1]:
                                #     file = open('ratio.txt','a')
                                #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                                if (linksStrength[1] < linksStrength[0] * thresholdRejected):  # then it means that the link does not exist
                                    se = segment._links[endOfSegment][n2]
                                    se.print_complete()
                                    segment._links[endOfSegment][n2].remove_end_of_link(segment._otherEndOfLinks[endOfSegment][n2], segment, endOfSegment)
                                    segment.remove_end_of_link(endOfSegment, segment._links[endOfSegment][n2], segment._otherEndOfLinks[endOfSegment][n2])
                                    se.print_complete()
                                    
                                elif (linksStrength[1] < linksStrength[0] * thresholdAccepted):  # then it's not clear, the link is freezed
                                    segment.freeze(endOfSegment)
    
                            else:
                                #     file = open('ratio.txt','a')
                                #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                                if linksStrength[0] < linksStrength[1] * thresholdRejected:  # then that the link does not exist
                                    segment._links[endOfSegment][n1].remove_end_of_link(segment._otherEndOfLinks[endOfSegment][n1], segment, endOfSegment)
                                    segment.remove_end_of_link(endOfSegment, segment._links[endOfSegment][n1], segment._otherEndOfLinks[endOfSegment][n1])
                                    
                                elif linksStrength[0] < linksStrength[1] * thresholdAccepted:  # then it's not clear, the link is freezed
                                    segment.freeze(endOfSegment)
                        n2+=1

    return listOfSegments

#merge_contigs looks at choices end of segments by endofsegment, and duplicates all the necessary contigs
def merge_contigs(listOfSegments, copiesnumber):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends) : once it is duplicated, it is locked

    #look at both ends of each segment sequentially
    for segment in listOfSegments:
        
        for endOfSegment in range(2) : 
            
            l = segment.links[endOfSegment]
            if len(l) > 1 and not segment.freezed[endOfSegment]:
                
                startMerging = not segment.locked
                for i in l:
                    startMerging = startMerging and (not i.locked)

                if startMerging:  # if nothing is locked for now
                    duplicate_around_this_end_of_contig(segment,endOfSegment,listOfSegments,copiesnumber)
                    
    #now that the duplicating is done, unlock all the segments
    print([i.locked for i in listOfSegments])
    print([i.names for i in listOfSegments])
    for segment in listOfSegments :
        segment.locked = False
    
    # now just merge all two contigs that are next to each other
    listOfSegments = merge_adjacent_contigs(listOfSegments)
    
    return listOfSegments, copiesnumber


def solve_ambiguities(listOfSegments, interactionMatrix, dist_law, stringenceReject, stringenceAccept, steps, copiesNumber = []):
        
    if copiesNumber == [] :
        copiesNumber = [1 for i in listOfSegments]

    listOfSegments = merge_adjacent_contigs(listOfSegments)

    for i in range(steps):

        print(str(i / steps * 100) + "% of solving ambiguities done\n")

        get_rid_of_bad_links(listOfSegments, interactionMatrix, dist_law, copiesNumber, stringenceReject, stringenceAccept)
            
        listOfSegments, copiesNumber = merge_contigs(listOfSegments, copiesNumber)
        
        #once all the contigs have been duplicated and merged, unfreeze everything so the cycle can start again
        for j in listOfSegments :
            j.unfreeze()
            j.print_complete()
        
        bf.export_to_GFA(listOfSegments, exportFile="tests/fake" + str(i) + ".gfa")
       # print('At step ', i, ' the energy is : ', score_output(listOfSuperContigs, links, [10000 for i in names], interactionMatrix))

    return listOfSegments


# links = bf.import_links('listsPython/links.csv')
# # print(links[1465])
# # infContigs = bf.read_info_contig('data/results/info_contigs.txt')
# interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')
# for i in range(len(interactionMatrix)):
#       interactionMatrix[i][i] = 0
# print('Loaded')

# # print(
# #     solve_ambiguities(
# #         [[], [4], [], [4], [1, 3], [6, 8], [5], [10], [5], [10], [7,9], []],
# #         [0, 1, 2, 3, 4, 5],
# #         [
# #             [1, 1, 1, 1, 0, 1],
# #             [1, 1, 1, 0, 1, 1],
# #             [1, 1, 1, 1, 1, 1],
# #             [1, 0, 1, 1, 1, 1],
# #             [0, 1, 1, 1, 1, 1],
# #             [1, 1, 1, 1, 1, 1]
# #         ]
# #     )
# # )

# newlinks, listOfSuperContigs, copiesnumber = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix, 2, 12)

# #print(copiesnumber)

# print('Now exporting')
# export_to_GFA(newlinks, listOfSuperContigs, links, copiesnumber, 'data/Assembly.fasta')

# print('Finished')

