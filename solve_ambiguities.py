#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import numpy as np
import basic_functions as bf
import scipy.integrate as integrate

from copy import deepcopy 

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
    
    commonContigs = []
    potentialCommonContigs = candidatesSegments[0].listOfContigs
    if supercontigsaretouching:
        if segment in candidatesSegments[0].links[0] :
            potentialCommonContigs += candidatesSegments[0].links[1]
        elif segment in candidatesSegments[0].links[1] :
            potentialCommonContigs += candidatesSegments[0].links[1]
        else :
            supercontigsaretouching = False
            print('WARNING in intensity_of_interactions : called like the contigs were touching but they were not')

    for contig in potentialCommonContigs:
        
        common = True

        for n, candidate in enumerate(candidatesSegments[1:]):

            allcontigs = candidate.listOfContigs.coopy()
            if supercontigsaretouching:
                
                if segment in candidate.links[0] :
                    allcontigs += candidate.links[1].copy()
                elif segment in candidate.links[1] :
                    allcontigs += candidate.links[0].copy()
                else :
                    supercontigsaretouching = False
                    print('WARNING in intensity_of_interactions : called like the contigs were touching but they were not')
     
            common = common and (contig in allcontigs)

        if common:
            commonContigs += [contig]

    # now we have a list of common contigs
    for candidate in candidatesSegments :
        if all(elem in commonContigs for elem in candidate.listOfContigs
        ):  # if all elements of sg are in commoncontigs, the algorithm cannot make a choice for now
            return [-1], [-1]

    #bestsignature decides what contig is most characterisic of the segment
    bestSignature = np.min([copiesnumber[x] for x in segment.listOfContigs])
    # we return for each supercontig its absolute score and its relative score (wihtout the common parts)
    absoluteScores = [0 for i in range(len(candidatesSegments))]
    relativeScores = [0 for i in range(len(candidatesSegments))]

    for c in candidatesSegments:
        
        absoluteScore, relativeScore, partial_area = c.interaction_with_contig(segment, interactionMatrix, dist_law, copiesnumber, commonContigs, bestSignature)

        absoluteScores.append(absoluteScore)
        relativeScores.append(relativeScore)
        
        if supercontigsaretouching:
            # print('partial area of contig ', candidatesSuperContigs[sg], ' : ', partial_area)
            relativeScores[-1] /= partial_area

    #print("*Common contigs : ", commonContigs)
    return absoluteScores, relativeScores


# here we look specifically at one contig and its immediate surroundings (can return -1 if fails in short loop)
def duplicate_around_this_end_of_contig(segment, endOfSegment, listOfSuperContigs, copiesnumber): #endOfSegment should be 0 if it's the left end and1 if it's the right end
    
    for i in segment.listOfContigs:
        copiesnumber[i] += len(segment.links[endOfSegment]) - 1

    # add all the new supercontigs
    for neighbor in segment.links[endOfSegment]:

        s.merge_two_segments(segment, endOfSegment, neighbor, listOfSuperContigs) #the merged segment is appended at the end of listOfSuperContigs

    # now delete the merged supercontigs
    # start by deleting the links that linked the merged supercontigs to the outside

    deletedContigs = [segment]
    otherEnd = 1 - endOfSegment

    for i, neighbor in enumerate(segment.links[otherEnd]):
    
        neighbor.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)
        

    for m, merged in enumerate(segment.links[endOfSegment]):
        
        if (len(merged.links[segment.otherEndOfLinks[endOfSegment][m]]) == 1):  # then the original copy is fully integrated in the supercontig
            deletedContigs += [merged]
            
            otherEnd = 1-segment.otherEndOfLinks[endOfSegment][m]
            for i, neighbor in enumerate(merged.links[otherEnd]):
                try:
                    neighbor.remove_end_of_link(merged.otherEndOfLinks[otherEnd][i], merged, otherEnd)
                except ValueError:  # that means we're in a small loop which whe can't solve
                    print("There is merging difficulty around the far end of "
                        + str(merged.namesOfContigs)+ " from "
                        + str(segment.namesOfContigs)+ " . Please check that there is indeed a loop there.")
                    return -1, -1, -1
                
        else:  # then the original contig still exists by itself, just delete the link going toward segment
            try:
                merged.remove_end_of_link(segment.otherEndOfLinks[endOfSegment][m], segment, endOfSegment)
            except ValueError:  # that means we're in a small loop which whe can't solve
                print("There is merging difficulty around the near end of "
                        + str(merged.namesOfContigs)+ " from "
                        + str(segment.namesOfContigs)+ " . Please check that there is indeed a loop there.")
                return -1, -1, -1

    #delete all segments that should be
    for i in deletedContigs:
        listOfSuperContigs.remove(i)

    return listOfSuperContigs, copiesnumber


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
        n.remove_end_of_link(neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor)

    # delete the ex-segments
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)

    return listOfSegments

#a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):
    
    n = len(listOfSegments)
    
    for se in range(n):
        segment = listOfSegments[se]
        for endOfSegment in range(2):
            if len(segment.links[endOfSegment]) == 1\
                and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]]) == 1:  # then merge
                if segment != segment.links[endOfSegment][0]:
                    listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)

    return listOfSegments

#get_rid_of_bad_links compare links using HiC contact informations when there is a choice and delete links that are not supported by HiC evidence
def get_rid_of_bad_links(listOfSegments,interactionMatrix,dist_law,copiesnumber,thresholdRejected,thresholdAccepted):

    # loop through all segments inspecting the robustness of all links.
    for segment in listOfSegments:
        
        endOfSegment = -1
        if len(segment.links[0]) > 1:
            endOfSegment = 0
        elif len(segment.links[1]) > 1 :
            endOfSegment = 1
            
        if endOfSegment >= 0 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
            
            mustBeDeleted = []
            mustBeDeletedOtherEnd = []

            # first : comparison pairwise of the links, those that should be deleted are stored in mustBeDeleted
            for n1 in range(len(segment.links[endOfSegment]) - 1):
                for n2 in range(n1 + 1, len(segment.links[endOfSegment])):
                    absoluteLinksStrength, linksStrength = intensity_of_interactions(segment,[listOfSegments[n1], listOfSegments[n2]]\
                                                                                     , listOfSegments, interactionMatrix, dist_law, copiesnumber, True)
                    if absoluteLinksStrength == [-1]: #means that the configuration does not enable the algorithm to compare the two interactions
                        segment.freezeNode(endOfSegment)
                        
                    else:
                        # print("I am ",listOfSuperContigs[int(endOfContig / 2)]," and here is the choice I have to make : ",
                        #     [listOfSuperContigs[int(links[endOfContig][neighbor1] / 2)],
                        #         listOfSuperContigs[int(links[endOfContig][neighbor2] / 2)]],linksStrength)

                        if linksStrength[0] > linksStrength[1]:
                            #     file = open('ratio.txt','a')
                            #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                            if (linksStrength[1] < linksStrength[0] * thresholdRejected):  # then that the link does not exist
                                mustBeDeleted.append(segment._links[endOfSegment][n2])
                                mustBeDeletedOtherEnd.append(segment._otherEndOfLinks[endOfSegment][n2])
                            elif (linksStrength[1] < linksStrength[0] * thresholdAccepted):  # then it's not clear, the link is freezed
                                segment.freeze(endOfSegment)

                        else:
                            #     file = open('ratio.txt','a')
                            #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                            if (
                                linksStrength[0] < linksStrength[1] * thresholdRejected
                            ):  # then that the link does not exist
                                mustBeDeleted.append(segment._links[endOfSegment][n1])
                                mustBeDeletedOtherEnd.append(segment._otherEndOfLinks[endOfSegment][n1])
                            elif (
                                linksStrength[0] < linksStrength[1] * thresholdAccepted
                            ):  # then it's not clear, the link is freezed
                                segment.freeze(endOfSegment)

            # second, the links that should be deleted are deleted
            for l in range(len(mustBeDeleted)) :
                segment.remove_end_of_link(endOfSegment, mustBeDeleted[l], mustBeDeletedOtherEnd[l])
                mustBeDeleted[l].remove_end_of_link(mustBeDeletedOtherEnd[l], segment, endOfSegment)

    return listOfSegments

#merge_contigs looks at choices end of segments by endofsegment, and duplicates all the necessary contigs
def merge_contigs(listOfSegments, copiesnumber):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends), so the other end is locked
    n = len(listOfSegments)

    #look at both ends of each segment sequentially
    for s in range(n):
        
        for endOfSegment in range(2) : 
            
            segment = listOfSegments[s]
            l = segment.links[endOfSegment]
            if len(l) > 1 and not segment.freezed[endOfSegment]:
                
                startMerging = not segment.locked
                for i in l:
                    startMerging = startMerging and (not i.locked)

                if startMerging:  # if nothing is locked for now
                    lsc, cn = duplicate_around_this_end_of_contig(segment,endOfSegment,deepcopy(listOfSegments),copiesnumber.copy())
                    
                    if lsc != -1:  # the function duplicate_... returns -1 when it fell on an unsolvable loop
    
                        listOfSegments, copiesnumber = deepcopy(lsc),deepcopy(cn)  # doing all these deepcopies takes time, we might want to reconsider to gain time

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
        
        # bf.export_to_GFA(
        #     links,
        #     listOfSuperContigs,
        #     copiesNumber,
        #     originalLinks,
        #     names,
        #     exportFile="tests/fake" + str(i) + ".gfa")
       # print('At step ', i, ' the energy is : ', score_output(listOfSuperContigs, links, [10000 for i in names], interactionMatrix))

    return listOfSegments, copiesNumber


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

