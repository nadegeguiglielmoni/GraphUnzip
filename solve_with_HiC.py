#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 18:06:07 2021

@author: rfaure

A file summarizing the algorithmic part of untangling the graph with HiC
"""

from determine_multiplicity import determine_multiplicity
from interaction_between_contigs import interactions_with_neighbors

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
def solve_with_HiC(segments, interactionMatrix, names, copiesnumber={}, confidentCoverage=True):
    
    if copiesnumber == {} :
        for segment in segments :
            for name in segment.names :
                copiesnumber[name] = 1
    
    #normalize the interactionMatrix
    normalInteractions = normalize(interactionMatrix)
    
    #print("Interactions of contig 147 : ", [normalInteractions[names["edge_147"], i] for i in range(len(names))])
    
    #determine all single-copy contigs
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links considering the topography of the graph
    refCoverage, multiplicities = determine_multiplicity(segments, names, supported_links)
    
    
    haploidContigs = []
    for se, s in enumerate(segments) :
    
        if multiplicities[se] == 1 :
        #to be deemed haploid, a segment must have at most one connection at each of its end plus be equally or less covered thant neighboring contigs
            links = s.links
            if round(s.depth/refCoverage) <= 1 :
                haploidContigs.append(s)
            # elif len(links[0]) == 1 and len(links[1]) == 1 and s.depth <= links[0][0].depth and s.depth <= links[1][0].depth : 
            #     haploidContigs.append(s)
            # elif len(links[0]) == 1 and len(links[1]) == 0 and s.depth <= links[0][0].depth : 
            #     haploidContigs.append(s)
            # elif len(links[0]) == 0 and len(links[1]) == 1 and s.depth <= links[1][0].depth : 
            #     haploidContigs.append(s)
        
    haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in haploidContigs :
        haploidContigsNames[s.full_name()] = index
        index += 1
                
    #determine explicitely all knots of the graph
    list_of_knots, list_of_neighbors, knotOfContig, haploidContigs, haploidContigsNames = determine_list_of_knots(segments, multiplicities, haploidContigs, haploidContigsNames, normalInteractions, names)
    
    #try to know what haploid contigs go together
    contacts = [[] for i in range(len(list_of_knots))]
    #while not sure :
    solvedKnots, rien1, rien2, sure = match_haploidContigs(segments, names, normalInteractions, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames, copiesnumber)
    
    untangled_paths = find_paths(contacts, segments, list_of_knots, solvedKnots, haploidContigsNames, haploidContigs, interactionMatrix, names, confidentCoverage)
    
    segments = untangle_knots(untangled_paths, segments)
    


#input: a graph, in the form of a list of segments and a list of haploid conigs
#output : a list of all independant knots of the graph. A new list of haploidContigs where each haploid contig actually has contacts with neighboring haploid contigs
def determine_list_of_knots(segments, multiplicities, haploidContigs, haploidContigsNames, interactionMatrix, names) :
    
    #first compute the list of neighbors for all haploidContigs, checking if each contig is informative or not
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    notInformative = set()
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames)
        
        #check if the contig is informative and if it does not loop on itself
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
            #     while True :
            #         r=0
            if all([i==0 for i in interactions]) :
                notInformative.add(end//2)
          
            #now check if it does not loop on itself
            for neighbor in listOfNeighbors[end] :
                if neighbor//2 == end//2 :
                    notInformative.add(end//2)
                    
            for neighbor in listOfNeighbors[end-1] :
                if neighbor//2 == end//2 :
                    notInformative.add(end//2)
            
    #remove all uninformative haploidContigs from haploidContigs
    #print("Here are all the uninformative contigs I should remove : ", [haploidContigs[i].names for i in notInformative])
    # while True :
    #     pass
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
                # print("I found the knot of ", s.full_name(), " : ", [haploidContigs[i//2].full_name() for i in listOfKnots[-1]])
    
    # for i in range(len(listOfNeighbors)) :
    #     print("Neighbors of contig ", haploidContigs[i//2].names, " are : ", [haploidContigs[j//2].names for j in listOfNeighbors[i]])
    
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
def match_haploidContigs(segments, names, interactionMatrix, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames, copiesnumber) :
    
    sure_haploids = True
    not_actually_haploid = [] #a list of contigs that have no contacts with neighbors among the haploidContigs (those contigs are useless)
    solvedKnots = []
    
    for k, knot in enumerate(list_of_knots) :
        
        if len(knot) > 1 :
            preferred_contact = {} #preferred_contact associates to an end all the ends that points toward it
            knot_solved = True
            
            for e, end in enumerate(knot) :
                
                index = haploidContigsNames[haploidContigs[end//2].full_name()]
                #print("Looking at interaction from contig ", haploidContigs[end//2].full_name())
                interactions = []
                for neighbor in list_of_neighbors[end] :
                    #print("Interaction with neighbor ", haploidContigs[neighbor//2].full_name())
                    
                    total = 0
                    for name1 in haploidContigs[end//2].names :
                        for name2 in haploidContigs[neighbor//2].names :
                            total += interactionMatrix[names[name1], names[name2]]
                    interactions += [ total ]

                
                interactions = interactions_with_neighbors(haploidContigs[end//2], end%2, [haploidContigs[i//2] for i in list_of_neighbors[end]], [i%2 for i in list_of_neighbors[end]], segments, interactionMatrix, names, copiesnumber)
                if interactions == [-1] :
                    print ("Did not manage to compute interactions from ", haploidContigs[end//2].names, " to ", [haploidContigs[i//2].names for i in list_of_neighbors[end]])
                
                print("Looking at interaction from contig ", haploidContigs[end//2].full_name(), " and here are its interactions: ", interactions)
                
                m = max(interactions)
                
                if m > 0 :
                    bestIdx = interactions.index(max(interactions))
                    
                    if list_of_neighbors[end][bestIdx] not in preferred_contact :
                        preferred_contact[list_of_neighbors[end][bestIdx]] = end  
                    # else : #it means two contigs point to the same contig, which is a problem
                    #     #knot_solved = False
                    #     #print("Contig ", haploidContigs[list_of_neighbors[end][bestIdx]//2].names, " does not look so haploid")
                    #     sure_haploids = False
                    
                    if end < list_of_neighbors[end][bestIdx] : #to put each link only once
                        contacts[k] += [(end, list_of_neighbors[end][bestIdx])] 
                    
                else :
                    print("I can't do anything about contig ", haploidContigs[end//2].names)
                    sure_haploids = False
                    not_actually_haploid += [end//2]
                    knot_solved = False

                    
            if knot_solved :
                solvedKnots += [k]
                print("We solved a knot:")
                for contact in contacts[k] :
                    print(haploidContigs[contact[0]//2].names, " -> ", haploidContigs[contact[1]//2].names)
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
            

    #print("Contacts : ", contacts)
            
    return solvedKnots, reliable_haploid_contigs, reliable_haploid_contigsNames, sure_haploids #return the updated list of haploid contigs

#a function normalizing an interaction matrix by iteratively dividing the rows and the columns by their norms. Additionally sets to 0 all values on the diagonal
#input : un-normalized matrix
#output : a new normalized matrix
def normalize(matrix) :
    
    matrixNow = matrix.tocoo()
    newMatrix = sparse.dok_matrix(matrix.shape)
    
    newMatrix = sparse.dok_matrix(matrix.shape)
    
    for steps in range(5) : #do 5 rounds of dividing
    
        sums = [0 for i in range(matrix.shape[0])]
        for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data):
            sums[r] += m
        for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data) :
            newMatrix[r,c] = m/sums[r]
          
        matrixNow = newMatrix.tocoo()
        
        sums = [0 for i in range(matrix.shape[0])]
        for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data):
            sums[c] += m
        for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data) :
            newMatrix[r,c] = m/sums[c]
            
        matrixNow = newMatrix.tocoo()
                
    #end by normalizing on the rows, because we will query on the rows
    sums = [0 for i in range(matrix.shape[0])]
    for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data):
        sums[r] += m
    for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data) :
        if r != c: #sets 0 on the diagonal
            newMatrix[r,c] = m/sums[r]
         
    # values = []
    # for r, c, m in zip(matrixNow.row, matrixNow.col, matrixNow.data):
    #     values += [m]
        
    # plt.hist(values, bins=50)
    # plt.xlim([0.05, 0.8])
    # plt.ylim([0,1000])
    # plt.show()
        
    print("Finished normalizing the interaction matrix")
    return newMatrix

#input : list of ends of haploid contigs we want to link, list of the solved knots we can work on (a partially resolved knot is no good)
#output : paths between the linked segments
def find_paths(contacts, segments, knots, solvedKnots, haploidContigsNames, haploidContigs, interactionMatrix, names, confidentCoverage) :
    
    untangled_paths = []
    
    for k in solvedKnots :
        
        untangled_paths += [[]]
        
        knot = knots[k]
        #the first step is to list all the ways to answer the contacts, i.e. for each path the set of decisions that can be made to link the two extremities
        #a set of decisions for 1 path comes in the form of a dictionnary of (contig, endOfContig) : listOfNeighbor. For example (contig1, 0) : [0,2] means that you can go from the extremity 0 of contig1 to its neighbor number 0 or 2 when coming from endOfPath2 and trying to go to endOfPath1
        alldecisions = [] #list of all sets of decisions for all paths of the knot
        touchedContigs = {} #a dict associating to each contig the list of paths it finds itself on
        for p, path in enumerate(contacts[k]) :
            
            decisions = find_decisions_on_path(haploidContigs[path[0]//2], path[0]%2, haploidContigs[path[1]//2], path[1]%2, haploidContigsNames, touchedContigs, p)
            alldecisions.append(decisions)
            
            # print("The knot is : ", [haploidContigs[i//2].names for i in knot], k)
            # print("Here are the decisions I can make to go from contig ", haploidContigs[path[0]//2].full_name(), " to ", haploidContigs[path[1]//2].full_name())
            # print(decisions)
    
        # print("touchedContigs : ", touchedContigs, "\n")
        
        repartitionOfContigs = dispatch_contigs(touchedContigs, contacts[k], interactionMatrix, names, haploidContigs, confidentCoverage)
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
        multiplicity = max([len(intercontig.links[0]), len(intercontig.links[1]), multiplicity])
        
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
            
        if 'edge_111' in intercontig.names :
                    print("Here is the repartition of contig 111 ", repartitionOfContigs[intercontig], " among the paths ", contacts)
                    
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
                        
                        newpath = best_path_coming_from_there[(segment, end)][1] + '<>'[1-otherEnd] + neighbor.full_name()
                        newdict = {i:best_path_coming_from_there[(segment,end)][2][i] for i in repartitionOfContigs.keys()}
                        newdict[neighbor] += 1
                        best_path_coming_from_there[(neighbor, 1-otherEnd)] = (potentialScore, newpath, newdict)
        
            list_of_new_paths_to_check = new_list_of_new_paths_to_check
        
        resultPaths += [best_path_coming_from_there[(segment1, 1-end1)][1]]
        
        print("Best path to go from ", segment2.names, " to ", segment1.names, " is ", resultPaths[-1])
    
    return resultPaths


def untangle_knots(untangled_paths, segments)    :
    
    fullnames = {}
    for s, seg in enumerate(segments) :
        fullnames[seg.full_name()] = s
        
    toDelete = set() #set of contigs to delete at the end
      
    for knot in range(len(untangled_paths)) :
        
        #compute how many times each intermediary contig is present in the final untangling
        numberofcopies = {}
        for path in untangled_paths[knot] :
            contigs = split('[><]' , path)
            del contigs[0]
            for contigName in contigs :
                if contigName in numberofcopies :
                    numberofcopies[contigName] += 1
                else :
                    numberofcopies[contigName] = 1
            
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
                    toDelete.add(fullnames[contig.full_name()])
                    
                elif c>0 and c == len(contigs) -1 : #then do not duplicate the segment, but do link it
                
                    end1 = '><'.index(orientations[c])
                    end0 = '<>'.index(orientations[c-1])
                    #print("The neighbors kept in mind to be neighbors of ", contig.names, " are ", [i.names for i in borderLinks[p][1]], ". I'm looking for ", segments[fullnames[contigs[c-1]]].names, " in path ", path)
                    idxNeighbor = borderLinks[p][1].index(segments[fullnames[contigs[c-1]]])
                    add_link(contig, end1, segments[newContigsIndices[c-1]], end0, borderCIGARs[p][1][idxNeighbor])
 
                    
    #delete all the contigs that were integrated
    newsegments = []
    for s, seg in enumerate(segments) :
        if s not in toDelete :
            newsegments += [seg]
        else :
            seg.cut_all_links()
    
    segments = newsegments
            
    return segments
            

        
    
    
    
    
    
    
    
