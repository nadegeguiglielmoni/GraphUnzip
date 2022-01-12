#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 18:06:07 2021

@author: rfaure

A file summarizing the algorithmic part of untangling the graph with HiC
"""

from determine_multiplicity import determine_multiplicity
from interaction_between_contigs import interactions_with_neighbors
from interaction_between_contigs import compute_commonContigs
from finish_untangling import merge_adjacent_contigs
from finish_untangling import break_up_chimeras

from segment import delete_link
from segment import Segment
from segment import add_link

from scipy import sparse
from scipy.sparse import isspmatrix_csr
from random import random
#import matplotlib.pyplot as plt
import numpy as np
from re import split #to split on several characters (<> here)
from re import findall  #to find all numbers in a mixed number/letters string (such as 31M1D4M)


#input : segments and the interactionMatrix of Hi-C contacts. Optionnaly, a list of haploid contigs obtained from the long reads algorithm.
def solve_with_HiC(segments, interactionMatrix, names, copiesnumber={}, confidentCoverage=True, verbose = False):
    
    if copiesnumber == {} :
        for segment in segments :
            for name in segment.names :
                copiesnumber[name] = 1
    
    #normalize the interactionMatrix
    normalInteractions = normalize(interactionMatrix, verbose)
    
    #print("Interactions of contig 147 : ", [normalInteractions[names["edge_147"], i] for i in range(len(names))])
    
    #determine the single-copy contigs that will serve as anchor for our algorithm
    haploidContigs = []
    
    ##compute the average depth of single-copy contigs
    totalDepth = 0
    totalLength = 0
    if confidentCoverage :
        for s in segments :
            if len(s.links[0]) <= 1 and len(s.links[1]) <= 1  : 
                totalDepth += s.depth * s.length
                totalLength += s.length
    refCoverage = totalDepth/totalLength            
    
    ##give a first estimation of the haploid contigs
    for se, s in enumerate(segments) :
    
        #to be deemed haploid, a segment must have at most one connection at each of its end plus be equally or less covered thant neighboring contigs
        links = s.links
        if round(s.depth/refCoverage) <= 1 and confidentCoverage:
            m1, m2 = 0, 0
            if len(links[0]) > 0 :
                m1 = max([i.depth for i in links[0]])
            if len(links[1]) > 0 :
                m2 = max([i.depth for i in links[1]])
            if s.depth < 1.5 * max(m1, m2) :
                haploidContigs.append(s)
        
        elif not confidentCoverage :
            if len(links[0]) <= 1 and len(links[1]) <= 1  : 
                haploidContigs.append(s)
    
    ##do not take as haploid anchors all contigs which are implicated in a knot by both ends
    haploidContigs_set = set(haploidContigs)
    toDelete = []
    for segment in haploidContigs_set :
        common_contigs = compute_commonContigs(segment, [segment, segment], [0,1], min(3*segment.length, 1000000))
        if len(common_contigs) > len(segment.names) :
            toDelete += [segment]
    for i in toDelete :
        haploidContigs_set.remove(i)
    haploidContigs = list(haploidContigs_set)
    
    
    #now all haploid contigs have been determined    
    haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in haploidContigs :
        haploidContigsNames[s.full_name()] = index
        index += 1
      
    solvedKnots = [0]  
    go_on = True
    limit = 5
    limit_counter = 0
    while go_on and limit_counter < limit: 
        
        limit_counter += 1
        
        #determine explicitely all knots of the graph
        list_of_knots, list_of_neighbors, knotOfContig, haploidContigs, haploidContigsNames = determine_list_of_knots(segments, haploidContigs, haploidContigsNames, normalInteractions, names, verbose)
        
        #try to know what haploid contigs go together
        contacts = [[] for i in range(len(list_of_knots))]
        #while not sure :
        
        solvedKnots, rien1, rien2, sure = match_haploidContigs(segments, names, normalInteractions, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames, copiesnumber, verbose)  
        untangled_paths = find_paths(contacts, segments, list_of_knots, solvedKnots, haploidContigsNames, haploidContigs, normalInteractions, names, confidentCoverage, verbose)
        
        segments, haploidContigs, haploidContigsNames, go_on = untangle_knots(untangled_paths, segments, haploidContigs)
        
    
    #segments = break_up_chimeras(segments, names, interactionMatrix, 1000000)
        
    
    return segments

#input: a graph, in the form of a list of segments and a list of haploid conigs
#output : a list of all independant knots of the graph. A new list of haploidContigs where each haploid contig actually has contacts with neighboring haploid contigs
def determine_list_of_knots(segments, haploidContigs, haploidContigsNames, interactionMatrix, names, verbose = False) :
    
    #first compute the list of neighbors for all haploidContigs, checking if each contig is informative or not
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    notInformative = set()
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = list(set(find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames))) #list(set()) to have all elements only once
        
        #check if the contig is informative
        if end % 2 == 1 :
            
            interactions = []
            for neighbor in listOfNeighbors[end] :
                total = 0
                for name1 in haploidContigs[end//2].names :
                    for name2 in haploidContigs[neighbor//2].names :
                        total += interactionMatrix[names[name1], names[name2]]
                interactions += [ total ]
            for neighbor in listOfNeighbors[end-1] :
                total = 0
                for name1 in haploidContigs[end//2].names :
                    for name2 in haploidContigs[neighbor//2].names :
                        total += interactionMatrix[names[name1], names[name2]]
                interactions += [ total ]
            
            # if 'edge_26' in haploidContigs[end//2].names :
            #     print("That's edge_26: ", interactions, " ", [haploidContigs[i//2].names for i in listOfNeighbors[end-1]],  [haploidContigs[i//2].names for i in listOfNeighbors[end]])

            if all([i==0 for i in interactions]) :
                notInformative.add(end//2)
          
            
    #remove all uninformative haploidContigs from haploidContigs
    #print("Here are all the uninformative contigs I should remove : ", [haploidContigs[i].names for i in notInformative])

    index = 0
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    
    for i in range(len(haploidContigs)) :
        
        if i not in notInformative :
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].full_name()] = index
            index += 1
            
    haploidContigs = reliable_haploid_contigs
    haploidContigsNames = reliable_haploid_contigsNames
        
    #recompute listOfNeighbors
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames)
        #print("Neighbors of ", haploidContigs[end//2].full_name(), " : ", [haploidContigs[i//2].full_name() for i in listOfNeighbors[end]])

    #Now move on to listing the knots
    
    listOfKnots = []
    knotOfContig = [-1 for i in range(len(haploidContigs)*2)]
    endOfSegmentAlreadySeen = [False for i in range(len(haploidContigs)*2)]
    
    bothEndsInTheSameKnot = set()
    
    for se, s in enumerate(haploidContigs):
        
        for end in range(2) :
            
            if not endOfSegmentAlreadySeen[se*2+end] : #means we're looking at a new knot
                
                knot = [se*2+end]
                
                for segIdx in knot : #knot will increase during the loop, it's normal
                    extension = listOfNeighbors[segIdx]
                    #print("Extensions of ", s.full_name(), " : ", [haploidContigs[i//2].full_name() for i in extension])
                    endOfSegmentAlreadySeen[segIdx] = True
                    knotOfContig[segIdx] = len(listOfKnots)

                    for ex in extension :
                        if not endOfSegmentAlreadySeen[ex] :
                            knot += [ex]
                            knotOfContig[ex] = len(listOfKnots)
                            endOfSegmentAlreadySeen[ex] = True
                            #print("Extending knot from ", s.full_name(), " with: ", haploidContigs[ex//2].full_name())
                            
                knot.sort() 
                removeThisPairOfEnds = []
                for e in range(1,len(knot)) :
                    if knot[e]//2 == knot[e-1]//2 :
                        bothEndsInTheSameKnot.add(knot[e]//2)
                        removeThisPairOfEnds.append(e)
                for toRemove in removeThisPairOfEnds[::-1] :
                    del knot[toRemove]
                    del knot[toRemove-1]
                listOfKnots += [knot]
                # print("I found the knot of ", s.full_name(), " : ", [haploidContigs[i//2].full_name() for i in listOfKnots[-1]])
    
    #once again recompute the reliable list of haploid contigs (those with both ends in the same knot eliminated)
    index = 0
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    
    for i in range(len(haploidContigs)) :
        
        if i not in bothEndsInTheSameKnot :
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].full_name()] = index
            index += 1
            
    haploidContigs = reliable_haploid_contigs
    haploidContigsNames = reliable_haploid_contigsNames
    
    #recompute one last time listOfneighbors
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames)
        
        #forbid l-loops. If it is the only thing possible from end, check end as non-informative (it's at least non-haploid)
        for n, neighbor in enumerate(listOfNeighbors[end]) :
            if neighbor//2 == end//2 and neighbor%2 == end%2 :
                del listOfNeighbors[end][n]
                if len(listOfNeighbors[end]) == 0 :
                    notInformative.add(end//2)
                    
    #recompute one last time the list of knots
    listOfKnots = []
    knotOfContig = [-1 for i in range(len(haploidContigs)*2)]
    endOfSegmentAlreadySeen = [False for i in range(len(haploidContigs)*2)]
    for se, s in enumerate(haploidContigs):
        
        for end in range(2) :
            
            if not endOfSegmentAlreadySeen[se*2+end] : #means we're looking at a new knot
                
                knot = [se*2+end]
                
                for segIdx in knot : #knot will increase during the loop, it's normal
                    extension = listOfNeighbors[segIdx]
                    #print("Extensions of ", s.full_name(), " : ", [haploidContigs[i//2].full_name() for i in extension])
                    endOfSegmentAlreadySeen[segIdx] = True
                    knotOfContig[segIdx] = len(listOfKnots)

                    for ex in extension :
                        if not endOfSegmentAlreadySeen[ex] :
                            knot += [ex]
                            knotOfContig[ex] = len(listOfKnots)
                            endOfSegmentAlreadySeen[ex] = True
                            #print("Extending knot from ", s.full_name(), " with: ", haploidContigs[ex//2].full_name())
                            
                listOfKnots += [knot]
    

    #print("List of knots: ", [[haploidContigs[i//2].names for i in j] for j in listOfKnots])

    
    return listOfKnots, listOfNeighbors, knotOfContig, haploidContigs, haploidContigsNames

#a recursive function to find to what haploid contig a contig can be linked
def find_neighbors(segment, end, segmentsAlreadyTraversed, haploidContigsNames) :
    
    res = []
    
    for n, neighbor in enumerate(segment.links[end]) :

        #print(neighbor.full_name())
        if (neighbor.full_name(), 1-segment.otherEndOfLinks[end][n]) not in segmentsAlreadyTraversed :
            
            segmentsAlreadyTraversed.add((neighbor.full_name(), 1-segment.otherEndOfLinks[end][n]) )
            
            if neighbor.full_name() in haploidContigsNames :
                res += [2*haploidContigsNames[neighbor.full_name()] + segment.otherEndOfLinks[end][n]]
            
            else:
                res += find_neighbors(neighbor, 1-segment.otherEndOfLinks[end][n], segmentsAlreadyTraversed, haploidContigsNames)
      
    return res

    
#input: a list of haploid contigs and their neighbors and  an interaction matrix
#output: contacts, a matrix matching pairwise the ends of the haploid contigs
def match_haploidContigs(segments, names, interactionMatrix, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames, copiesnumber, verbose) :
    
    sure_haploids = True
    not_actually_haploid = [] #a list of contigs that have no contacts with neighbors among the haploidContigs (those contigs are useless)
    solvedKnots = []
    
    for k, knot in enumerate(list_of_knots) :
        
        if len(knot) > 1  :
            preferred_contact = {} #preferred_contact associates to an end all the ends that points toward it
            knot_solved = True
            
            for e, end in enumerate(knot) :
                
                index = haploidContigsNames[haploidContigs[end//2].full_name()]

                interactions = interactions_with_neighbors(haploidContigs[end//2], end%2, [haploidContigs[i//2] for i in list_of_neighbors[end]], [i%2 for i in list_of_neighbors[end]], segments, interactionMatrix, names, copiesnumber, verbose = verbose)
                if interactions == [-1] and verbose :
                    print ("Did not manage to compute interactions from ", haploidContigs[end//2].names, " to ", [haploidContigs[i//2].names for i in list_of_neighbors[end]])
                
                #if 'edge_128' in haploidContigs[end//2].names :
                    
                #print("Looking at interaction from contig ", haploidContigs[end//2].full_name(), " and here are its interactions: ", interactions, k)
                
                m = max(interactions)
                
                if m > 0 :
                    bestIdx = interactions.index(max(interactions))
                    
                    if list_of_neighbors[end][bestIdx] not in preferred_contact :
                        preferred_contact[list_of_neighbors[end][bestIdx]] = {end}
                    else :
                        preferred_contact[list_of_neighbors[end][bestIdx]].add(end)

                    contacts[k] += [(min(end, list_of_neighbors[end][bestIdx]), max(end, list_of_neighbors[end][bestIdx]))] 
                    
                else :
                    # print("I can't do anything about contig ", haploidContigs[end//2].names)
                    sure_haploids = False
                    not_actually_haploid += [end//2]
                    knot_solved = False

                    
            if knot_solved :
                
                #remove the links present twice
                contacts[k] = list(set(contacts[k]))
                
                #sometimes big contigs tend to be linked spuriously with each other : delete all contacts that are not absolutely necessary
                for contact in contacts[k][::-1] :
                    if contact[0] in preferred_contact and len(preferred_contact[contact[0]]) > 1 and contact[1] in preferred_contact and len(preferred_contact[contact[1]]) > 1 :
                        contacts[k].remove(contact)
                
                solvedKnots += [k]
                if verbose :
                    print("We solved knot ", [haploidContigs[i//2].names for i in knot])
                    for contact in contacts[k] :
                        print(haploidContigs[contact[0]//2].names, " -> ", haploidContigs[contact[1]//2].names)
            if verbose :
                print()
         
    #now get rid of non-haploid and uninformative contigs
    
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    not_actually_haploid = list(set(not_actually_haploid))
    not_actually_haploid.sort()
    
    index = 0
    indexNot = 0
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(not_actually_haploid) and i == not_actually_haploid[indexNot] :
            indexNot += 1
        else :
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
            

    if verbose :
        print("Finished matching haploid contigs, now we'll move on to the paths linking them")
            
    return solvedKnots, reliable_haploid_contigs, reliable_haploid_contigsNames, sure_haploids #return the updated list of haploid contigs

#a function normalizing an interaction matrix by iteratively dividing the rows and the columns by their norms. Additionally sets to 0 all values on the diagonal
#input : un-normalized matrix
#output : a new normalized matrix
def normalize(matrix, verbose) :
    
    if verbose :
        print("Normalizing the interaction matrix")
    W = matrix.tocsr()
    W = W.astype('float64')
    
    for rounds in range(10) :
        if verbose :
            print("Round ", rounds, "/10", end='\r')
        for i in range(W.shape[0]):
                row_sum = W.data[W.indptr[i]:W.indptr[i+1]].sum()
                if row_sum != 0:
                    W.data[W.indptr[i]:W.indptr[i+1]] /= row_sum
             
        W = W.transpose()
        W = W.tocsr()
        for i in range(W.shape[0]):
                row_sum = W.data[W.indptr[i]:W.indptr[i+1]].sum()
                if row_sum != 0:
                    W.data[W.indptr[i]:W.indptr[i+1]] /= row_sum
        W = W.transpose()
        
    #normalize one last time, as the interactions will be queried by row
    for i in range(W.shape[0]):
        row_sum = W.data[W.indptr[i]:W.indptr[i+1]].sum()
        if row_sum != 0:
            W.data[W.indptr[i]:W.indptr[i+1]] /= row_sum
            
    if verbose :
        print("Finished normalizing the interaction matrix")
        
    return W.todok()

#input : list of ends of haploid contigs we want to link, list of the solved knots we can work on (a partially resolved knot is no good)
#output : paths between the linked segments
def find_paths(contacts, segments, knots, solvedKnots, haploidContigsNames, haploidContigs, interactionMatrix, names, confidentCoverage, verbose = False) :
    
    untangled_paths = []
    
    for k in solvedKnots :
        
        untangled_paths += [[]]
        
        knot = knots[k]
        #the first step is to list all the ways to answer the contacts, i.e. for each path the set of decisions that can be made to link the two extremities
        #a set of decisions for 1 path comes in the form of a dictionnary of (contig, endOfContig) : listOfNeighbor. For example (contig1, 0) : [0,2] means that you can go from the extremity 0 of contig1 to its neighbor number 0 or 2 when coming from endOfPath2 and trying to go to endOfPath1
        alldecisions = [] #list of all sets of decisions for all paths of the knot
        touchedContigs = {} #a dict associating to each contig the list of paths it finds itself on
        if verbose :
                print("Expliciting all the paths for knot ", k)
        for p, path in enumerate(contacts[k]) :
            
            decisions = find_decisions_on_path(haploidContigs[path[0]//2], path[0]%2, haploidContigs[path[1]//2], path[1]%2, haploidContigsNames, touchedContigs, p)
            alldecisions.append(decisions)
            
            # print("The knot is : ", [haploidContigs[i//2].names for i in knot], k)
            # print("Here are the decisions I can make to go from contig ", haploidContigs[path[0]//2].full_name(), " to ", haploidContigs[path[1]//2].full_name())
            # print(decisions)
    
        # print("touchedContigs : ", touchedContigs, "\n")
        if verbose :
                print("Dispatching all the intermediary contigs")
        repartitionOfContigs = dispatch_contigs(touchedContigs, contacts[k], interactionMatrix, names, haploidContigs, confidentCoverage)
        if verbose :
            print("Finding dynamically the best paths to satisfy all constraints")
        untangled_paths[-1] = find_best_paths(alldecisions, repartitionOfContigs, segments, contacts[k], haploidContigs)
        
        #untangled_paths[-1] += find_best_paths(alldecisions, interactionMatrix, names, segments, contacts[k], haploidContigs, confidentCoverage)
        
        # print("Here are the untangled paths I got : ", untangled_paths)
        # print()
    
    return untangled_paths

#input: two ends of a path
#output : the list of all possible decisions we can make to go from segment1 to segment2, as well as a completed touchedContigs, i.e. with indexOfPath added on all possible contigs reached by this path
def find_decisions_on_path(segment1, end1, segment2, end2, haploidContigsNames, touchedContigs, indexOfPath) :
    
    touchingSegment1 = set() #a set of tuple (contig ID, end of contig) from which it is possible to reach segment 1 
    reachable_from_segment1(segment1, end1, haploidContigsNames, touchingSegment1)
        
    decisions = {} #decisions is a dict : for each end of each intermediary contig, there is a list of possible (neighbor,end), gives you all the neighbor indices you can choose to go from 2 to 1
    backtrack_from_segment2(segment2, end2, segment1, end1, touchingSegment1, decisions, touchedContigs, indexOfPath)
      
    return decisions

#input : a knot and an end of segment
#output : all the segments that can be reached from segment1 and from which end (in form of a set of tuple (ID, end))
def reachable_from_segment1(segment1, end1, haploidContigsNames, reachable) :
        
    for n, neighbor in enumerate(segment1.links[end1]) :
        
        otherEnd = segment1.otherEndOfLinks[end1][n]
        
        if (neighbor.full_name(), otherEnd) not in reachable and neighbor.full_name() not in haploidContigsNames :
        
            reachable.add((neighbor.full_name(), otherEnd))
            reachable_from_segment1(neighbor, 1-otherEnd, haploidContigsNames, reachable)
        
#input : all the contigs reachable from segment1
#output : all the decisions with which you can go from segment2 to segment 1, and touchedContigsCompleted with indexOfPath
def backtrack_from_segment2(segment2, end2, segment1, end1, touchingSegment1, decisions, touchedContigs, indexOfPath) :
    
    decisions[(segment2, end2)] = []
    
    for n, neighbor in enumerate(segment2.links[end2]) :
        
        otherEnd = segment2.otherEndOfLinks[end2][n]
        
        
        if (neighbor.full_name(), 1-otherEnd) in touchingSegment1 :
            
            decisions[(segment2, end2)] += [n]
            if neighbor not in touchedContigs :
                touchedContigs[neighbor] = {indexOfPath}
            else :
                touchedContigs[neighbor].add(indexOfPath)
            
            if (neighbor, 1-otherEnd) not in decisions :
                backtrack_from_segment2(neighbor, 1-otherEnd, segment1, end1, touchingSegment1, decisions, touchedContigs, indexOfPath)
                
        if neighbor.ID == segment1.ID and otherEnd == end1 :
            decisions[(segment2, end2)] += [n]
  
#input : a list of contacts and intermediary contigs
#output : a repartition of all intermediary contigs between the different paths
def dispatch_contigs(touchedContigs, contacts, interactionMatrix, names, haploidContigs, confidentCoverage) :

    repartitionOfContigs= {} #a dict associating each contig to a dictionnary giving the repartition of this contig
    
    #first deal with intermediary contigs
    refCoverage = np.mean([haploidContigs[i[0]//2].depth for i in contacts]+[haploidContigs[i[1]//2].depth for i in contacts]) #refCoverage is the average of the coverage of the haploid contigs bordering the knot

    for intercontig in touchedContigs.keys() :

        repartitionOfContigs[intercontig] = {}        

        #how many times this contig should appear ? :
        multiplicity = 0
        if confidentCoverage :
            multiplicity = round(intercontig.depth / refCoverage)
        multiplicity = max([ min(len(intercontig.links[0]), len(intercontig.links[1]), 2), multiplicity])
        
        #compute the interaction of this contig with each path it can be on
        interaction_with_path = [0 for i in range(len(contacts))]
        for p in touchedContigs[intercontig] :
            for subcontig in intercontig.names :
                for subco in haploidContigs[contacts[p][0]//2].names :
                    interaction_with_path[p] += interactionMatrix[names[subcontig], names[subco]]
                for subco in haploidContigs[contacts[p][1]//2].names :
                    interaction_with_path[p] += interactionMatrix[names[subcontig], names[subco]]
             
        #give an estimate of how it should be spread between the different paths
        totalInteraction = np.sum(interaction_with_path)

        for p in touchedContigs[intercontig] :
            if totalInteraction > 0 :
                estimate = round(multiplicity*interaction_with_path[p]/totalInteraction) #the number of times we expect to find this intermediary in this path
            else :
                estimate = -1
                
            repartitionOfContigs[intercontig][p] = estimate
            
        # if 'edge_111' in intercontig.names :
        #             print("Here is the repartition of contig 111 ", repartitionOfContigs[intercontig], " among the paths ", contacts)
                    
    #then dispatch the border contig (1 in each path in which they belong)
    for p, path in enumerate(contacts) :
        
        for extremity in range(2) :
            segment = haploidContigs[path[extremity]//2]
            end = path[extremity]%2
        
            if segment not in repartitionOfContigs :
                repartitionOfContigs[segment] = {}
            
            repartitionOfContigs[segment][p] = 1
                    
    return repartitionOfContigs
  
#input : a list of decisions that can be made for each path and the ideal repartition of intermediary contigs
#output : a proposition of a set of good paths to untangle the knot 
def find_best_paths(alldecisions, repartitionOfContigs, segments, contacts, haploidContigs) :
    
    resultPaths = []
    for p, path in enumerate(contacts) :
        
        segment1 = haploidContigs[path[0]//2]
        end1 = path[0]%2
        
        segment2 = haploidContigs[path[1]//2]
        end2 = path[1]%2
        
        best_path_coming_from_there = {} #a dictionary with the best score on each intermediary contig, the path it came from and the number of occurences of each contig on this path
        for intermediary in repartitionOfContigs.keys() :
            best_path_coming_from_there[(intermediary,0)] = (-10000, '', {i:0 for i in repartitionOfContigs.keys()})
            best_path_coming_from_there[(intermediary,1)] = (-10000, '', {i:0 for i in repartitionOfContigs.keys()})
        best_path_coming_from_there[(segment1, 1-end1)] = (-10000, '', {i:0 for i in repartitionOfContigs.keys()})
        best_path_coming_from_there[(segment2, end2)] = (0, '<>'[end2]+segment2.full_name(), {i:0 for i in repartitionOfContigs.keys()}) #this is the start point of the exploration
            
        list_of_new_paths_to_check = [(segment2, end2)] #start from there
        while len(list_of_new_paths_to_check) > 0 :
            
            new_list_of_new_paths_to_check = []
            list_of_check_scores = [] # a list of the scores of the paths to check, because we want to set a limit to the number of paths we explore, since this number grows combinatorily in very complex regions
            for (segment, end) in list_of_new_paths_to_check :
                
                # print("Here is the end, to go from ", segment2.names, " to ", segment1.names, " : ", segment.names, " ", end)
                # print("All the keys of decisions : ",[(i[0].names, i[1]) for i in alldecisions[p].keys()])
                for n in alldecisions[p][(segment,end)] :
                    
                    neighbor = segment.links[end][n]
                    otherEnd = segment.otherEndOfLinks[end][n]
                    
                    if repartitionOfContigs[neighbor][p] != -1 : #a score of -1 meant a contig we did not manage to dispatch
                        good = repartitionOfContigs[neighbor][p] - best_path_coming_from_there[(segment,end)][2][neighbor] - 0.5 #positive if we lack this contig in this path, negative elsewhise
                        potentialScore = best_path_coming_from_there[(segment,end)][0] + good/abs(good) 
                    else : #if we have no information on this contig, the contig is weightless
                        potentialScore = best_path_coming_from_there[(segment,end)][0]
                    
                    if potentialScore > best_path_coming_from_there[(neighbor, 1-otherEnd)][0] :
                        
                        if neighbor.ID != segment1.ID or otherEnd != end1 :
                            new_list_of_new_paths_to_check += [(neighbor, 1-otherEnd)]
                            list_of_check_scores += [potentialScore]
                        
                        newpath = best_path_coming_from_there[(segment, end)][1] + '<>'[1-otherEnd] + neighbor.full_name()
                        #newdict = {i:best_path_coming_from_there[(segment,end)][2][i] for i in repartitionOfContigs.keys()}
                        newdict = best_path_coming_from_there[(segment,end)][2].copy()
                        newdict[neighbor] += 1
                        best_path_coming_from_there[(neighbor, 1-otherEnd)] = (potentialScore, newpath, newdict)
        
            best_paths = list(range(len(list_of_check_scores)))
            best_paths.sort(reverse = True, key = lambda x : list_of_check_scores[x])
            list_of_new_paths_to_check = [new_list_of_new_paths_to_check[i] for i in best_paths[:1000]] #only explore the 1000 best paths, it is largely sufficient
        
        resultPaths += [best_path_coming_from_there[(segment1, 1-end1)][1]]
        
        #print("Best path to go from ", segment2.names, " to ", segment1.names, " is ", resultPaths[-1])
    
    return resultPaths

#input : a list of paths for each knot
#output : untangled graph in the updated segments
def untangle_knots(untangled_paths, segments, haploidContigs)    :
    
    fullnames = {}
    for s, seg in enumerate(segments) :
        fullnames[seg.full_name()] = s
        
    numberOfSegments_start = len(segments)
        
    toDelete = set() #set of contigs to delete at the end
    endSolved = [[False, False] for i in range(len(segments))] #an array stocking True if the end of this segment is at the end of an haploidcontig on a solved path
    go_on = False #a boolean switching to True if anything is modified
      
    for knot in range(len(untangled_paths)) :
        
        #compute how many times each intermediary contig is present in the final untangling
        numberofcopies = {}
        for path in untangled_paths[knot] :
            contigs = split('[><]' , path)
            del contigs[0]
            orientations = "".join(findall("[<>]", path))
            for contigName in contigs :
                if contigName in numberofcopies :
                    numberofcopies[contigName] += 1
                else :
                    numberofcopies[contigName] = 1
                    
            endSolved[fullnames[contigs[0]]]['<>'.index(orientations[0])] = True
            endSolved[fullnames[contigs[-1]]]['><'.index(orientations[0])] = True
            
        #delete all links of the haploid contigs pointing toward the knot (we will reconstruct the paths later). However, carefully keep the CIGARs and what they point to, we don't want to lose them
        borderCIGARs = []
        borderLinks = []
        
        ##first save all the CIGARs
        for p, path in enumerate(untangled_paths[knot]) :
            
            contigs = split('[><]' , path)
            orientations = "".join(findall("[<>]", path))
            del contigs[0]
            
            borderCIGARs += [[[],[]]]
            borderLinks += [[[],[]]]
            for extremity in range(2) :
                segment = segments[fullnames[contigs[-extremity]]]
                end = '<>'.index(orientations[0])
                if extremity == 1 :
                    end = '><'.index(orientations[-1])
                
                for n in range(len(segment.links[end])): #delete all the links right of the contig, the only good one will be reestablished later (actually, only delete the links to haploidContigs)
                    borderCIGARs[-1][extremity] += [segment.CIGARs[end][n]]
                    borderLinks[-1][extremity] += [segment.links[end][n]]
        
        ##then delete the links
        for p, path in enumerate(untangled_paths[knot]) :
            
            contigs = split('[><]' , path)
            orientations = "".join(findall("[<>]", path))
            del contigs[0]
            
            for extremity in range(2) :
                segment = segments[fullnames[contigs[-extremity]]]
                end = '<>'.index(orientations[0])
                if extremity == 1 :
                    end = '><'.index(orientations[-1])
                
                while len(segment.links[end]) > 0 : #delete all the links right of the contig, the only good one will be reestablished later (actually, only delete the links to haploidContigs)
                    neighbor = segment.links[end][0]
                    delete_link(segment, end, neighbor, segment.otherEndOfLinks[end][0])
        
        #now reconstruct the path of contigs leading from one haploid contig to the other side of the knot
        for p, path in enumerate(untangled_paths[knot]) :
                                    
            #the first and last contig of path are haploid, the rest are intermediary contigs
            contigs = split('[><]' , path)
            orientations = "".join(findall("[<>]", path))
            del contigs[0] #because it is always ''
            
            newContigsIndices = [fullnames[contigs[0]]]
            
            for c, contigName in enumerate(contigs) :
                
                contig = segments[fullnames[contigName]]
                
                
                if c > 0 and c < len(contigs)-1 :
                    
                    go_on = True
                    
                    newSegment = Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i/numberofcopies[contigName] for i in contig.depths])
                    segments.append(newSegment)
                    newContigsIndices += [len(segments) - 1]
                    
                    #add the link to form the new bridge
                    end1 = '><'.index(orientations[c])
                    end0 = '<>'.index(orientations[c-1])
                    if c == 1 :
                        idxNeighbor = borderLinks[p][0].index(segments[fullnames[contigs[1]]])
                        CIGAR = borderCIGARs[p][0][idxNeighbor]
                    else :
                        idxNeighbor = contig.links[end1].index(segments[fullnames[contigs[c-1]]]) 
                        CIGAR= contig.CIGARs[end1][idxNeighbor]
                        
                    add_link(segments[-1], end1, segments[newContigsIndices[c-1]], end0, CIGAR)
                    
                    #now that we have created another contig, the old one should be deleted at the end
                    toDelete.add(segments[fullnames[contig.full_name()]])
                    
                elif c>0 and c == len(contigs) -1 : #then do not duplicate the segment, but do link it
                
                    end1 = '><'.index(orientations[c])
                    end0 = '<>'.index(orientations[c-1])
                    #print("The neighbors kept in mind to be neighbors of ", contig.names, " are ", [i.names for i in borderLinks[p][1]], ". I'm looking for ", segments[fullnames[contigs[c-1]]].names, " in path ", path)
                    idxNeighbor = borderLinks[p][1].index(segments[fullnames[contigs[c-1]]])
                    add_link(contig, end1, segments[newContigsIndices[c-1]], end0, borderCIGARs[p][1][idxNeighbor])
 
            
    #now look if there aren't supposedly haploid segments that are linked twice or more at their ends
    for s in range(numberOfSegments_start) :
              
        contig = segments[s]
        for end in range(2) :
            
            #if it is solved at only one end, duplicate the contig
            #print(endSolved[s])
            # if 'edge_120' in contig.names :
            #     print("Is it solved : ", end, endSolved[s][end], len(contig.links[end]) , )
            if endSolved[s][end] and len(contig.links[end]) > 1 and len(contig.links[1-end]) == 0 and contig.depth > 1.5*np.sum([i.depth for i in contig.links[end]]):
                
                contig.divide_depths(len(contig.links[end]))
                
                for n, neighbor in enumerate(contig.links[end]):
                    if n > 0 :
                        newSegment = Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i for i in contig.depths])
                        add_link(newSegment, end, neighbor, contig.otherEndOfLinks[end][n], contig.CIGARs[end][n])
                        delete_link(contig, end, neighbor, contig.otherEndOfLinks[end][n])
                        
                        #add also the links from the other side
                        for n2, neighbor2 in enumerate(contig.links[1-end]) :
                            add_link(newSegment, 1-end, neighbor2, contig.otherEndOfLinks[1-end][n2], contig.CIGARs[1-end][n2])
                        
                        segments.append(newSegment)
 
    #delete all the contigs that were integrated
    newsegments = []
    for s, seg in enumerate(segments) :
        if seg not in toDelete :
            newsegments += [seg]
        else :
            seg.cut_all_links()
    
    newsegments = merge_adjacent_contigs(newsegments)
    
    #now go through all contigs and create a new haploidContigs list
    haps = set(haploidContigs)
    pastSegments = set([i.ID for i in segments])

    stillHaploids = [] 
    for s, seg in enumerate(newsegments) :
        if (seg in haps or seg not in pastSegments) and (len(seg.links[0]) == 1 or len(seg.links[1]) == 1 or seg.length > min(np.min([i.length for i in seg.links[0]], initial=seg.length),np.min([i.length for i in seg.links[1]], initial=seg.length))) : #if seg is not in pastsegment, it has been created thus is haploid (if it is long enough or has 1 link at one of its ends)
            
            stillHaploids += [seg]
            
    #now all haploid contigs have been determined 
    stillHaploids.sort(key= lambda x: x.length, reverse = True)
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in stillHaploids :
        haploidContigsNames[s.full_name()] = index
        index += 1
            
    return newsegments, stillHaploids, haploidContigsNames, go_on
            

        
    
    
    
    
    
    
    