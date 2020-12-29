#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import numpy as np
import input_output as io
from bisect import bisect_left #to look through sorted lists

from copy import deepcopy
import os

from transform_gfa import check_segments
import segment as s
from segment import Segment

#this function detects and breaks up long (>length) chimeric contigs
def break_up_chimeras(segments, names, interactionMatrix, length) :
    
    allsegments = []
    allXs = []
    for s, segment in enumerate(segments) :
        
        if segment.length > length :
            
            interactions = []
            X = []
            
            for axis in range(1, len(segment.names)) :
                
                interaction = 0
                
                for nameLeft in segment.names[axis:] :
                    for nameRight in segment.names[:axis] :
                        
                        nameLeft = nameLeft.split('-')[0]
                        nameRight = nameRight.split('-')[0]
                        
                        interaction += interactionMatrix[names[nameLeft], names[nameRight]]
                        
                interactions += [interaction]
                
                if axis > 1 :
                    X += [ X[-1] + segment.lengths[axis-1]/segment.length ]
                else :
                    X = [segment.lengths[0]/segment.length]
                
            allsegments += [interactions]
            
            #plt.plot (X, interactions)
            
            inSlump = False
            localMinimums = []
            for axis in range(1, len(interactions)-1) :
                
                if interactions[axis] < 0.7*np.max(interactions[:axis]) and interactions[axis] < 0.7*np.max(interactions[axis:]):
                    if not inSlump :
                        inSlump = True
                        
                    localMinimums += [axis]
                        
                
                else :
                    if inSlump :
                        inSlump = False
                        
                        loin = [interactions[i] for i in localMinimums].index( np.min([interactions[i] for i in localMinimums]) )
                        
                        print('Breaking up contig ', segment.names, ' between ', segment.names[localMinimums[loin]-1], ' and ', segment.names[localMinimums[loin]], ' because it looks like a chimeric contig')
                                
                        #Now break the contig where it should
                        newSegment1, newSegment2 = segment.break_contig(localMinimums[loin])
                        segments[s] = newSegment1
                        segments.append(newSegment2)
                        
                        localMinimums = []

        
    #plt.show()    
    return segments

# this function measures the intensity of interactions between one supercontig
# and several candidate, including without taking account of the common parts
# of the supercontigs
# It also weighs the interaction with the length of a supercontig, so that a
# very long candidate spercontig is not seen as having a lot of connexion just
# because it is long
def intensity_of_interactions(
    segment,
    candidatesSegments,
    listOfTouchingEnds,
    listOfSegments,
    interactionMatrix,
    names, 
    copiesnumber,
    depthOfCommonContigs = 2, #meaning by default you take neighbors of neighbors in common contigs. (value 0,1,2)
    debugDir = '',
):
    for candidate in candidatesSegments :
        if candidate == segment : #small loop, don't solve that !
            return [-1], [-1], True #the True value does not matter here
    
    ##first compute all contigs common to all candidates, to take them out
    if depthOfCommonContigs > 0 :
            commonContigs, neighborsOfNeighborsUsed = compute_commonContigs(segment, candidatesSegments, listOfTouchingEnds, depthOfCommonContigs)
            if debugDir != '':
                f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
                f.write('common contigs: ' + str(commonContigs)+'\n')
                f.close()
    

    ##Now compute the score of each candidates    

    # bestsignature decides what contig is most characterisic of the segment
    bestSignature = np.min([copiesnumber[x] for x in segment.names])
    # return for each supercontig its absolute score and its relative score (wihtout the common parts)
    
    depthFound = True
    depth = 2
    while depthFound == True : #do the process again without neighborOfneighbor if neighbors of neighbors are only used sometimes but not always
        absoluteScores = []
        relativeScores = []
        returnRelativeScore = True
        depthFound = False

        for c in candidatesSegments:
    
            if depth == 2 :
                    
                absoluteScore, relativeScore, depthHere = c.interaction_with_contigs(segment, interactionMatrix, names, copiesnumber, commonContigs, bestSignature, neighborsOfNeighborsUsed)
        
                if depthHere == 1 and depth == 2 :
                    depth = 1
                    depthFound = True
                    
            else :
                absoluteScore, relativeScore, depthHere = c.interaction_with_contigs(segment, interactionMatrix, names, copiesnumber, commonContigs, bestSignature, False)
                
            if all([i in commonContigs for i in c.names]) :
                returnRelativeScore = False #if a contig is entirely contained in commonContigs, don't say that his score is 0, that would lead to errors
    
            absoluteScores.append(absoluteScore)
            relativeScores.append(relativeScore)
            
    if all([i==0 for i in relativeScores]) :
        returnRelativeScore = False # if all elements of candidate are in commoncontigs, relative intensity cannot be determined, you have to do with absolute intensity
    
    # if 'edge_229' in segment.names:    
        # print('At contig ', segment.names, ' choosing between ',  [i.names for i in candidatesSegments], ' and the result is ', relativeScores, absoluteScores)
        # #print('Best signature : ', bestSignature, ' and the signatures are : ', [copiesnumber[x] for x in segment.names])
        # print('Common contigs : ', commonContigs, '\n')
    
    if returnRelativeScore :
        return absoluteScores, relativeScores, neighborsOfNeighborsUsed
    else :
        return absoluteScores, [-1], neighborsOfNeighborsUsed

def compute_commonContigs(segment, candidatesSegments, listOfTouchingEnds, depth) :

    commonContigs = []
    if depth == 0 :
        return commonContigs
        
    if depth == 2 :
    #first compute the list of common contigs counting neighbors and neighbors of neighbors of segment
        potentialCommonContigs = candidatesSegments[0].names.copy()
    
        for i in candidatesSegments[0].links[1-listOfTouchingEnds[0]] :
            potentialCommonContigs += i.names
    
    
        for contig in potentialCommonContigs:
            
    
            for n in range(1, len(candidatesSegments)):
    
                candidate = candidatesSegments[n]
                presentInNeighbor = contig in candidate.names
    
                for i in candidate.links[1-listOfTouchingEnds[n]] :                      
                    presentInNeighbor = presentInNeighbor or (contig in i.names)
                    
    
            if presentInNeighbor:
                commonContigs += [contig]
                
        #check if that list of common contigs is not too big
        neighborOfNeighborActivated = True
        for c in candidatesSegments :
            if all([elem in commonContigs for elem in c.names]) : #this is not good, commoncontigs is too big
                neighborOfNeighborActivated = False
        
        if neighborOfNeighborActivated :
            return commonContigs, True #the True value is to signifie that neighbors of neighbors were used
    
    
    
    #recompute common contigs but only with neighbors
    commonContigs = []
    potentialCommonContigs = candidatesSegments[0].names.copy()
        
    for contig in potentialCommonContigs:
    
        common = True

        for n, candidate in enumerate(candidatesSegments[1:]):

            common = common and (contig in candidate.names)

        if common:
            commonContigs += [contig]
        
    return commonContigs, False #the False value is to signifie that neighbors of neighbors were not used


# small function to look in a sorted list l if x is present in logarithmic time
def isPresent(l, x):
    i = bisect_left(l, x)
    if i != len(l) and l[i] == x:
        return True
    return False


# here we look specifically at one contig and its immediate surroundings (can
# return -1 if fails in short loop)
def duplicate_around_this_end_of_contig(
    segment, endOfSegment, listOfSuperContigs, copiesnumber
):  # endOfSegment should be 0 if it's the left end and1 if it's the right end

    if (
        segment in segment.links[endOfSegment]
    ):  # if a segment loops on itself, another module would be needed
        return 0
    if any(
        segment.links[endOfSegment].count(i) >= 2 for i in segment.links[endOfSegment]
    ):  # if segment.links[endOfSegment] has two copies of the same segment, it means one link going towards each end, that is not solvable
        return 0

    for i in segment.names:
        copiesnumber[i] += len(segment.links[endOfSegment]) - 1

    # add all the new supercontigs
    for neighbor in segment.links[endOfSegment]:
        s.merge_two_segments(
            segment, endOfSegment, neighbor, listOfSuperContigs
        )  # the merged segment is appended at the end of listOfSuperContigs

    # now delete the merged supercontigs
    # start by deleting the links that linked the merged supercontigs to the outside

    deletedContigs = [segment.ID]
    otherEnd = 1 - endOfSegment

    for i, neighbor in enumerate(segment.links[otherEnd]):

        neighbor.remove_end_of_link(
            segment.otherEndOfLinks[otherEnd][i], segment, otherEnd
        )

    for m, merged in enumerate(segment.links[endOfSegment]):

        if (
            len(merged.links[segment.otherEndOfLinks[endOfSegment][m]]) == 1
        ):  # then the original copy is fully integrated in the supercontig
            deletedContigs.append(merged.ID)

            otherEnd = 1 - segment.otherEndOfLinks[endOfSegment][m]
            for i, neighbor in enumerate(merged.links[otherEnd]):
                try:
                    neighbor.remove_end_of_link(
                        merged.otherEndOfLinks[otherEnd][i], merged, otherEnd
                    )
                except ValueError:  # that means we're in a small loop which we can't solve
                    print(
                        "There is merging difficulty around the far end of "
                        + str(merged.names)
                        + " from "
                        + str(segment.names)
                        + " . Please check that there is indeed a loop there."
                    )
                    return 0

        else:  # then the original contig still exists by itself, just delete the link going toward segment
            try:
                merged.remove_end_of_link(
                    segment.otherEndOfLinks[endOfSegment][m], segment, endOfSegment
                )
            except ValueError:  # that means we're in a small loop which whe can't solve
                print(
                    "There is merging difficulty around the near end of "
                    + str(merged.names)
                    + " from "
                    + str(segment.names)
                    + " . Please check that there is indeed a loop there."
                )
                return 0

    # delete all segments that should be
    deletedContigs.sort()
    for i in range(len(listOfSuperContigs) - 1, -1, -1):
        h = listOfSuperContigs[i].ID
        if isPresent(
            deletedContigs, h
        ):  # in other words, if h is in deletedContigs (written like that because it has logarithmic efficiency)
            del listOfSuperContigs[i]

    # lock all the segments that have been duplicated, so that they are not duplicated by both ends
    segment.lockNode(endOfSegment)

# similar to the function above, but simpler: put in one supercontig two smaller supercontig linked by a link unambinguous at both ends
def merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments):

    if len(segment.links[endOfSegment]) != 1:
        print("ERROR : trying to merge simply two contigs that cannot be merged simply")
        return -1, -1

    neighbor = segment.links[endOfSegment][0]
    endOfSegmentNeighbor = segment.otherEndOfLinks[endOfSegment][0]

    if len(neighbor.links[endOfSegmentNeighbor]) != 1 :
        print('ERROR : trying to merge simply two contigs that cannot be merged simply')
        return -1,-1
        
    if neighbor == segment :  # then do not merge a contig with itself
        return -1, -1

    # add the new segment
    s.merge_two_segments(segment, endOfSegment, neighbor, listOfSegments)

    # delete links going towards the two ex-segments
    otherEnd = 1 - endOfSegment
    otherEndNeighbor = 1 - endOfSegmentNeighbor

    for i, n in enumerate(segment.links[otherEnd]) :
        n.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)
        
    for i, n in enumerate(neighbor.links[otherEndNeighbor]) :
        #print('Removing ', neighbor.names, ' from ', n.names, ' and adding the new contig',listOfSegments[-1].names, ' at end ', neighbor.otherEndOfLinks[otherEndNeighbor][i])
        n.remove_end_of_link(neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor, otherEndNeighbor)

    # delete the ex-segments    
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)
    
    return listOfSegments

# a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):

    goOn = True
    while goOn:
        goOn = False
        for segment in listOfSegments:

            alreadyDidThisOne = False # if the segment is deleted when looking at its first end, you don't want it to look at its other end, since it does not exist anymore
            for endOfSegment in range(2):
                
                if not alreadyDidThisOne:
                    
                    if len(segment.links[endOfSegment]) == 1and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]])== 1:  # then merge
                        alreadyDidThisOne = True
                        if segment != segment.links[endOfSegment][0]:
                            goOn = True
                            listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)

    return listOfSegments

#merge_contigs looks at choices endofsegment by endofsegment, and duplicates all the necessary contigs
def merge_contigs(listOfSegments, copiesnumber, verbose = False):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends) : once it is duplicated, it gets locked

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
                    
                    if verbose :
                        print('Duplicating contig ', segment.names, ' around its end touching ', [i.names for i in listOfSegments])
                    
    #now that the duplicating is done, unlock all the segments
    for segment in listOfSegments :
        segment.locked = False
        
    # now just merge all two contigs that are next to each other
    listOfSegments = merge_adjacent_contigs(listOfSegments)
    
    return listOfSegments, copiesnumber
             
#function to freeze and sometimes solve small o-loops (works only if long reads are there)
def solve_small_loops(listOfSegments, names, repeats, lr_links, check_links) :
    
    for se in range(len(listOfSegments)) :
        
        segment = listOfSegments[se]
        if segment in segment.links[0] : #this is a small loop of length 0
            
            if segment.links[0].count(segment) == 1 : #this is a o-loop, let's flatten it
                
                if repeats != [] :
                    replications = 0
                    for contig in segment.names :
                        replications = max(replications, repeats[names[contig]])
                        
                    segment.flatten(replications)
                    #print('In solve_small_loops, flattening ', segment.names, segment.insideCIGARs)
                
                else :
                    segment.freeze(0)
                    segment.freeze(1)
        

        toRemove = []
        for n, neighbor in enumerate(segment.links[0]) : #trying to detect o-loops of length 1
            endOfLink = segment.otherEndOfLinks[0][n]
            index = s.find_this_link(neighbor, 1-endOfLink, segment.links[1], segment.otherEndOfLinks[1], warning = False) #returns -1 if it does not find anything
            
            if index != -1 : # then this is a o-loop of length 1
                    
                    if lr_links != [] :
                        cA0 = segment.names[0]
                        oA0 = (0 == segment.orientations[0])
                        cB0 = neighbor.names[-endOfLink]
                        oB0 = (neighbor.orientations[-endOfLink] == endOfLink)
                        
                        if not (cA0, oA0, cB0, oB0) in lr_links and not (cB0, oB0, cA0, oA0) in lr_links :
                            toRemove += [(segment, 0, neighbor, endOfLink)]
                        
                        endOfLink2 = segment.otherEndOfLinks[1][index]
                        cA1 = segment.names[-1]
                        oA1 = (1 == segment.orientations[-1])
                        cB1 = neighbor.names[-endOfLink2]
                        oB1 = (neighbor.orientations[-endOfLink2] == endOfLink2)
                        if not (cA1, oA1, cB1, oB1) in lr_links and not (cB1, oB1, cA1, oA1) in lr_links:
                            toRemove += [(segment, 1, neighbor, endOfLink2)]
                        
                        
        for i in toRemove :
            # print('In o-loops : removing link from ', i[0].names, i[1], 'to ', i[2].names, i[3])
            # print('Links from ', i[0].names, i[1], ' : ', [j.names for j in i[0].links[i[1]]], i[0].otherEndOfLinks[i[1]])
            # print('Links from ', i[2].names, i[3], ' : ', [j.names for j in i[2].links[i[3]]], i[2].otherEndOfLinks[i[3]])
            i[0].remove_end_of_link(i[1], i[2], i[3])
            i[2].remove_end_of_link(i[3], i[0], i[1])
                
                
            
def solve_l_loops(segments, lr_links): #l-loops occur when one end of a contig is in contact with both end of another contig or when a contig is linked to itself at one end
    
    for segment in segments :
        for endOfSegment in range(2) :
            toRemove = []
                    
            for n in range(len(segment.links[endOfSegment])-1) :
                    
                if segment.links[endOfSegment][n].ID == segment.links[endOfSegment][n+1].ID and  segment.otherEndOfLinks[endOfSegment][n] == segment.otherEndOfLinks[endOfSegment][n+1]: #the two links going toward the same contig are next to each other because links are sorted
                    
                    #here we have a l-loop
                    neighbor = segment.links[endOfSegment][n]
                    lenToRemove = len(toRemove)
                    
                    if lr_links != [] :
                        #let's check if the two links of the l-loop are confirmed by long reads
                        
                        if neighbor.ID == segment.ID : #here we have a l-loop of length 0
                            cA = segment.names[-endOfSegment]
                            oA = (endOfSegment == segment.orientations[-endOfSegment])
                            if not (cA, oA, cA, not oA) in lr_links :
                                toRemove += [(segment, endOfSegment, segment, endOfSegment)]
                            
                        else : #l-loop of length 1
                        
                            cA0 = segment.names[-endOfSegment]
                            oA0 = (endOfSegment == segment.orientations[-endOfSegment])
                            cB0 = neighbor.names[-segment.otherEndOfLinks[endOfSegment][n]]
                            oB0 = (neighbor.orientations[-segment.otherEndOfLinks[endOfSegment][n]] == segment.otherEndOfLinks[endOfSegment][n])
                            if not (cA0, oA0, cB0, oB0) in lr_links and not (cB0, not oB0, cA0, not oA0) in lr_links :
                                toRemove += [(segment, endOfSegment, neighbor, int(oB0))]
                                    
                            cA1 = segment.names[-endOfSegment]
                            oA1 = (endOfSegment == segment.orientations[-endOfSegment])
                            cB1 = neighbor.names[-segment.otherEndOfLinks[endOfSegment][n+1]]
                            oB1 = (neighbor.orientations[-segment.otherEndOfLinks[endOfSegment][n+1]] == segment.otherEndOfLinks[endOfSegment][n+1])
                            if not (cA1, oA1, cB1, oB1) in lr_links and not (cB1, not oB1, cA1, not oA1) in lr_links:
                                toRemove += [(segment, endOfSegment, neighbor, int(oB1))]
                    
                    if len(toRemove) == lenToRemove : #means that no links could be taken out
                        segment.freeze(endOfSegment)
                        
            for i in toRemove :
                i[0].remove_end_of_link(i[1], i[2], i[3])
                i[2].remove_end_of_link(i[3], i[0], i[1])
        
#input : a graph. This function takes out all links that are not confirmed by long reads
def check_all_links(segments, lr_links) :
    
    
    for segment in segments :
        
        for endOfSegment in range(2) :
            
            toRemove = set()
            
            for n, neighbor in enumerate(segment.links[endOfSegment]) :
                
                
                endOfNeighbor = segment.otherEndOfLinks[endOfSegment][n]
                
                #define the link in the format of the lr_links
                link = (segment.names[-endOfSegment], (segment.orientations[-endOfSegment] == endOfSegment), neighbor.names[-endOfNeighbor], (neighbor.orientations[-endOfNeighbor] == endOfNeighbor))
                linkb = (link[2], link[3], link[0], link[1])
                
                                
                if not link in lr_links and not linkb in lr_links : #then the link is not confirmed by long reads
                    
                    toRemove.add((segment, endOfSegment, neighbor, endOfNeighbor))
            
            for i in toRemove :
                
                i[0].remove_end_of_link(i[1], i[2], i[3])
                i[2].remove_end_of_link(i[3], i[0], i[1])
    
#get_rid_of_bad_links compare links using HiC contact informations when there is a choice and delete links that are not supported by HiC evidence
def get_rid_of_bad_links(listOfSegments, interactionMatrix, lrInteractionMatrix, names, copiesnumber,thresholdRejected,thresholdAccepted, lr_links, debugDir = '', neighborsOfNeighbors = True, verbose = False, exhaustive = True):

    HiCmatrix = (interactionMatrix.count_nonzero() > 0) #a boolean value to tell if there is need to use the Hi-C interaction matrix

    #loop through all segments inspecting the robustness of all links.
    c = 0
    
    #then compute the intensity of interactions knowing the common contigs
    for segment in listOfSegments:
        
        c += 1

        for endOfSegment in range(2):
                
            if len(segment.links[endOfSegment]) >= 2 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
                
                    
                # comparison pairwise of the links, those that should be deleted are deleted
                    toRemove = []
                    for n1 in range(len(segment.links[endOfSegment]) - 1):
                        n2 = n1 + 1
                        while n2 < len(segment.links[endOfSegment]):
                            
                            d = 2
                            if not neighborsOfNeighbors :
                                d = 1
                            
                            #first compute using long reads    
                            absoluteLinksStrength, linksStrength, neighborsOfNeighborsUsed = intensity_of_interactions(segment, [segment.links[endOfSegment][n1], segment.links[endOfSegment][n2]],\
                                                                                             [segment.otherEndOfLinks[endOfSegment][n1], segment.otherEndOfLinks[endOfSegment][n2]],\
                                                                                             listOfSegments, lrInteractionMatrix, names, copiesnumber, depthOfCommonContigs = d, debugDir = debugDir)
                                
                            tmpls = linksStrength.copy()    
                                          
                            #if it is not enough, use Hi-C
                                    
                            if not exhaustive or ( (linksStrength == [-1] or (all([i>1 for i in linksStrength]) or all([i<=1 for i in linksStrength]))) and HiCmatrix ):
                                absoluteLinksStrength, linksStrength, neighborsOfNeighborsUsed = intensity_of_interactions(segment, [segment.links[endOfSegment][n1], segment.links[endOfSegment][n2]],\
                                                                                                [segment.otherEndOfLinks[endOfSegment][n1], segment.otherEndOfLinks[endOfSegment][n2]],\
                                                                                                listOfSegments, interactionMatrix, names, copiesnumber, depthOfCommonContigs = d, debugDir = debugDir)
                                
                                if not exhaustive and linksStrength != [-1] and tmpls != [-1]:
                                    linksStrength = [linksStrength[i] +tmpls[i] for i in range(len(linksStrength))]
                            
                            if debugDir != '' :
                                f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
                                f.write('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength)+ '\t'+str(absoluteLinksStrength)+'\n')
                                f.close()
                                
                            if verbose :
                                print('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength)+ '\t'+str(absoluteLinksStrength))
                            # if 'edge_229' in segment.names : 
                            
                            #     print('At 229, choosing between ', segment.links[endOfSegment][n1].names, segment.links[endOfSegment][n2].names, ' with these values : ', linksStrength, absoluteLinksStrength, neighborsOfNeighborsUsed)
                                
                            if not neighborsOfNeighborsUsed : #means that there are a lot of common contigs, a sort of knot
                                segment.freeze(endOfSegment)
                                #print('get_rid_of_bad_links, ...  freeznoding : ' + '\t'.join( ['_'.join(segment.links[endOfSegment][n1].names), '_'.join(segment.links[endOfSegment][n2].names)])+'\n')
                            
                            if linksStrength == [-1]: #means that the configuration does not enable the algorithm to compare the two interactions
                                segment.freezeNode(endOfSegment)                         

                            elif any([i>1 for i in linksStrength]): #the condition is to prevent too much duplicating if there is no mapping or almost  
                                
                                # if '262' in segment.names :
                                #     print('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength)+'\n')
                                #     print([i.names for i in segment.links[0]])
                                if linksStrength[0] > linksStrength[1]:
                                    if (linksStrength[1] <= linksStrength[0] * thresholdRejected) or (linksStrength[1] == 1 and linksStrength[0] > 2):  # then it means that the link does not exist
                                        if verbose :
                                            print('\nRemoving link from ', segment.links[endOfSegment][n2].names, ' to ', segment.names, '\n')
                                        if n2 not in toRemove :
                                            toRemove += [n2]
                                        
                                    elif (linksStrength[1] < linksStrength[0] * thresholdAccepted):  # then it's not clear, the link is freezed
                                        segment.freezeNode(endOfSegment)

                                else:
                                    if linksStrength[0] < linksStrength[1] * thresholdRejected or (linksStrength[0] == 1 and linksStrength[1] > 2):  # then decide that the link does not exist
                                        if verbose :
                                            print('\nRemoving link from ', segment.links[endOfSegment][n1].names, ' to ', segment.names, '\n')
                                        if n1 not in toRemove :
                                            toRemove += [n1]
                                        
                                        
                                    elif linksStrength[0] < linksStrength[1] * thresholdAccepted:  # then it's not clear, the link is freezed
                                        segment.freezeNode(endOfSegment)
                            else : #linksStrength <= [1,1]
                                segment.freezeNode(endOfSegment)
                                # print('get_rid_of_bad_links, ...  freeznoding2 : ' + '\t'.join( ['_'.join(segment.links[endOfSegment][n1].names), '_'.join(segment.links[endOfSegment][n2].names)])+'\n')

                            n2+=1
                            
                    #Remove all links that have been marked as removable
                    toRemove.sort()
                    toRemove.reverse()
                    
                    for n in toRemove :
                        segment._links[endOfSegment][n].remove_end_of_link(segment._otherEndOfLinks[endOfSegment][n], segment, endOfSegment)
                        segment.remove_end_of_link(endOfSegment, segment._links[endOfSegment][n], segment._otherEndOfLinks[endOfSegment][n])
                        
        
    return listOfSegments                        

def solve_ambiguities(listOfSegments, interactionMatrix, lrInteractionMatrix, names, stringenceReject, stringenceAccept, steps, copiesNumber = {}, repeats = [], lr_links = [], useNeighborOfNeighbor = True, debugDir = '', check_links = False, verbose = False):
        
    if debugDir != '' :
        if not os.isdir(debugDir) :
            os.mkdir(debugDir)
        f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'w')
        f.close()
    
    if copiesNumber == {} :
        for segment in listOfSegments :
            copiesNumber['_'.join(segment.names)] = 1
    
    if check_links :
        check_all_links(listOfSegments, lr_links) # check if all links there are present in the long reads and delete those who are not
    
    
    listOfSegments = merge_adjacent_contigs(listOfSegments)
    print('Merged adjacent contigs for the first time')            

   # s.check_if_all_links_are_sorted(listOfSegments)
    
    for i in range(steps):
        
        get_rid_of_bad_links(listOfSegments, interactionMatrix, lrInteractionMatrix, names, copiesNumber, stringenceReject, stringenceAccept,  lr_links, debugDir = debugDir, neighborsOfNeighbors = useNeighborOfNeighbor, verbose = verbose, exhaustive = check_links)
        
        solve_small_loops(listOfSegments, names, repeats, lr_links, check_links)
        
        solve_l_loops(listOfSegments, lr_links)
            
        print('Got rid of bad links')

        # for se in listOfSegments :
        #     if 'edge_357' in se.names :
        #         print ('Here is two : ', se.names, [i.names for i in se.links[0]], [i.names for i in se.links[1]], '\n')

        listOfSegments, copiesNumber = merge_contigs(listOfSegments, copiesNumber, verbose = verbose)
    

        #print('end of merge_contigs : ', [i.names for i in listOfSegments[names['262']].links[0]])

        # once all the contigs have been duplicated and merged, unfreeze everything so the cycle can start again
        for j in listOfSegments:
            j.unfreeze()
        
        print(str((i+1) / steps * 100) + "% of solving ambiguities done")
        
        if debugDir != '' :
            io.export_to_GFA(listOfSegments, exportFile = debugDir.strip('/')+'/'+'debug_gfa_step'+str(i)+'.gfa')
            f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
            f.write('Finished step '+ str(i)+ ' \n\n\n')
            f.close()
            
    #finish by breaking up long chimeras that can form sometimes
    
    if interactionMatrix.count_nonzero() > 0 :
        listOfSegments =  break_up_chimeras(listOfSegments, names, interactionMatrix, 100000)
        
    return listOfSegments, copiesNumber #return copiesNumber in case you want to run solve_ambiguities several times in a row



