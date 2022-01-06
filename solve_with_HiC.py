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
    
    untangle_knots(untangled_paths, segments)


#input: a graph, in the form of a list of segments and a list of haploid conigs
#output : a list of all independant knots of the graph. A new list of haploidContigs where each haploid contig actually has contacts with neighboring haploid contigs
def determine_list_of_knots(segments, multiplicities, haploidContigs, haploidContigsNames, interactionMatrix, names) :
    
    #first compute the list of neighbors for all haploidContigs, checking if each contig is informative or not
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    notInformative = []
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames)
        
        #check if the contig is informative and it it does not loop on itself
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
                notInformative += [end//2]
          
            #now check if it does not loop on itself
            for neighbor in listOfNeighbors[end] :
                if neighbor//2 == end//2 :
                    notInformative += [end//2]
                    
            for neighbor in listOfNeighbors[end-1] :
                if neighbor//2 == end//2 :
                    notInformative += [end//2]
            
    #remove all uninformative haploidContigs from haploidContigs
    print("Here are all the uninformative contigs I should remove : ", [haploidContigs[i].names for i in notInformative])
    # while True :
    #     pass
    index = 0
    indexNot = 0
    
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(notInformative) and i == notInformative[indexNot] :
            indexNot += 1
        else :
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
        for path in contacts[k] :
            
            decisions = find_decisions_on_path(haploidContigs[path[0]//2], path[0]%2, haploidContigs[path[1]//2], path[1]%2, haploidContigsNames)
            alldecisions.append(decisions)
            
            print("The knot is : ", [haploidContigs[i//2].names for i in knot], k)
            print("Here are the decisions I can make to go from contig ", haploidContigs[path[0]//2].full_name(), " to ", haploidContigs[path[1]//2].full_name(), " : ")
            print(decisions)
    
        untangled_paths[-1] += find_best_paths(alldecisions, interactionMatrix, names, segments, contacts[k], haploidContigs, confidentCoverage)
        
        print("Here are the untangled paths I got : ", untangled_paths)
        print()
    
    return untangled_paths

#input: two ends of a path
#output : the list of all possible decisions we can make to go from segment1 to segment2
def find_decisions_on_path(segment1, end1, segment2, end2, haploidContigsNames) :
    
    touchingSegment1 = set() #a set of tuple (contig ID, end of contig) from which it is possible to reach segment 1 
    reachable_from_segment1(segment1, end1, haploidContigsNames, touchingSegment1)
        
    decisions = {} #decisions is a dict : for each end of each intermediary contig, there is a list of possible (neighbor,end), gives you all the neighbor indices you can choose to go from 2 to 1
    backtrack_from_segment2(segment2, end2, segment1, end1, touchingSegment1, decisions)
      
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
#output : all the decisions with which you can go from segment2 to segment 1
def backtrack_from_segment2(segment2, end2, segment1, end1, touchingSegment1, decisions) :
    
    decisions[(segment2.full_name(), end2)] = []
    
    for n, neighbor in enumerate(segment2.links[end2]) :
        
        otherEnd = segment2.otherEndOfLinks[end2][n]
        
        
        if (neighbor.full_name(), 1-otherEnd) in touchingSegment1 :
            
            decisions[(segment2.full_name(), end2)] += [n]
            
            if (neighbor.full_name(), 1-otherEnd) not in decisions :
                backtrack_from_segment2(neighbor, 1-otherEnd, segment1, end1, touchingSegment1, decisions)
                
        if neighbor.ID == segment1.ID and otherEnd == end1 :
            decisions[(segment2.full_name(), end2)] += [n]
            
#input : a list of decisions that can be made for each path
#output : a proposition of a set of good paths to untangle the knot
def find_best_paths(alldecisions, interactionMatrix, names, segments, contacts, haploidContigs, confidentCoverage) :
    
    #use a gradient descent to find the best set of paths, starting with a naive algorithm thathas equiprobability for all decisions
    
    rules_of_decisions = [] #this will be a list of dictionnary : to each decision will be attached a probability to do this or that
    
    #initiate rules_of_decisions to equiprobability
    for decisions in alldecisions :
        
        probabilities = {}
        for choices in decisions.keys() :
            
            probabilities[choices] = [1/len(decisions[choices]) for i in range(len(decisions[choices]))]
            
        rules_of_decisions += [probabilities]
            
    # loop and change rules_of_decisions until we don't manage to do better anymore
    
    time_since_last_improval = 0
    bestScore = 0
    bestPaths = []
    number_of_tries = 10 #we launch 10 times the same rule of decision before making it evolve
    number_of_generation = 0
    while time_since_last_improval < 2 : #if you don't manage to imrpove the paths twice in a row, call it a day
    
        #initiate the variable that will keep in mind all the decisions
        decisions_made = []
        for decisions in alldecisions :
        
            score = {}
            for choices in decisions.keys() :
                
                score[choices] = [[] for i in range(number_of_tries)]
                
            decisions_made += [score]
            
        #initiate the variable that keeps in mind the score of each try
        scores = [0 for i in range(number_of_tries)]
        
        #now try different things and see how it works
        for attempt in range(number_of_tries) :
            
            paths = try_a_untangling(attempt, scores, decisions_made, alldecisions, rules_of_decisions, segments, contacts, haploidContigs, interactionMatrix, names, confidentCoverage) #paths there is in format [>contig1<contig456<contig3 , <contig23>contig67>contig888]
            if scores[attempt] > bestScore :
                bestPaths = paths
                bestScore = scores[attempt]
                time_since_last_improval = -1
        
        #modify the rule of decision according to what has worked
        for p, path in enumerate(rules_of_decisions) :
        
            for choices in path.keys() :
                
                if len(path[choices]) > 1 :
                
                    applicable_scores = [scores[i] for i in range(len(scores)) if decisions_made[p][choices][i] != []]
                    mean_score =  np.mean(applicable_scores)
                    
                    for attempt in range(number_of_tries) :
                        if decisions_made[p][choices][attempt] != [] :
                            
                            if scores[attempt] > mean_score :
                                for ch in decisions_made[p][choices][attempt] :
                                    path[choices][ ch ] += 1/number_of_tries * 1/np.sqrt(1+number_of_generation) / len(decisions_made[p][choices][attempt])
                                if path[choices][ decisions_made[p][choices][attempt] ] > 0.9 :
                                    path[choices][ decisions_made[p][choices][attempt] ] = 0.9
                            elif scores[attempt] < mean_score :
                                for ch in decisions_made[p][choices][attempt] :
                                    path[choices][ ch ] -= 1/number_of_tries * 1/np.sqrt(1+number_of_generation) / len(decisions_made[p][choices][attempt])
                                if path[choices][ decisions_made[p][choices][attempt] ] < 0.1 :
                                    path[choices][ decisions_made[p][choices][attempt] ] = 0.1
                    
                    #we've tweaked the probabilities, but we must adjust them so they sum to one :
                    sumProbas = np.sum(path[choices])
                    path[choices] = [i/sumProbas for i in path[choices]]
                    
        #now rule_of_decision has changed, re-try other paths in the next iteration
                
        number_of_generation += 1
        time_since_last_improval += 1
        
    return bestPaths

#input : a rule of decision and a knot (implicitely described by contacts and haploidContigs)
#output : a scored set of paths obtained with the rule of decision
def try_a_untangling(attempt, scores, decisions_made, alldecisions, rules_of_decisions, segments, contacts, haploidContigs, interactionMatrix, names, confidentCoverage) :
    
    paths = []
    intermediaryContigs = {} #dictionnary associating to each intermediary contig the number of times it has been used in each path
    
    for p, path in enumerate(contacts) :
        
        segment1 = haploidContigs[path[0]//2]
        end1 = path[0]%2
        
        segment2 = haploidContigs[path[1]//2]
        end2 = path[1]%2
        
        #print("Tentatively building a path from ", segment1.full_name(), " to ", segment2.full_name())
        
        paths += ['<>'[end2]+segment2.full_name()]
        
        #start from the end of the path and go back to the beginning, following alldecisions
        while  segment2.ID != segment1.ID or end2 != 1-end1 :
            
            #choose where to go next
            #print("I'm going from ", haploidContigs[path[1]//2].names, " to ", haploidContigs[path[0]//2].names, ", in ", segment2.names, " now, here are the decisions I can take : ", alldecisions[p][(segment2.full_name(), end2)], " p: ", p, (segment2.full_name(), end2))
            decision = np.random.choice(alldecisions[p][(segment2.full_name(), end2)], p = rules_of_decisions[p][(segment2.full_name(), end2)])
            decisions_made[p][(segment2.full_name(), end2)][attempt] += [decision]
                        
            newsegment2 = segment2.links[end2][decision]
            newend2 = 1-segment2.otherEndOfLinks[end2][decision]
            
            #print("Making the path from ", segment2.full_name(), " to ", newsegment2.full_name(), ", thanks to decision ", decision)
            
            segment2 = newsegment2
            end2 = newend2
            
            if segment2.ID in intermediaryContigs :
                if p in intermediaryContigs[segment2] :
                    intermediaryContigs[segment2][p] += 1
                else :
                    intermediaryContigs[segment2][p] = 1
            else :
                intermediaryContigs[segment2] = {p:1}
            
            paths[-1] += '<>'[end2]+segment2.full_name()
     
    scores[attempt] = score_this_attempt(paths, interactionMatrix, names, contacts, haploidContigs, intermediaryContigs, confidentCoverage)
    return paths
            
#input : a set of paths, an interaction matix, an inventory of where intermediary contigs are found and in how many times
#output : a score for this set of paths
def score_this_attempt(paths, interactionMatrix, names, contacts, haploidContigs, intermediaryContigs, confidentCoverage) :
    
    score = 0
    refCoverage = np.mean([haploidContigs[i[0]//2].depth for i in contacts]+[haploidContigs[i[1]//2].depth for i in contacts]) #refCoverage is the average of the coverage of the haploid contigs bordering the knot
    
    for intercontig in intermediaryContigs.keys() :
        
        interaction_with_path = [0 for i in range(len(contacts))]
        for p in range(len(contacts)) :
            
            #compute the interaction of this contig with each path
            for subcontig in intercontig.names :
                for subco in haploidContigs[contacts[p][0]//2].names :
                    interaction_with_path[p] += interactionMatrix[names[subcontig], names[subco]]
                for subco in haploidContigs[contacts[p][1]//2].names :
                    interaction_with_path[p] += interactionMatrix[names[subcontig], names[subco]]
            
        #how many times this contig should appear ? :
        multiplicity = 0
        if confidentCoverage :
            multiplicity = round(intercontig.depth / refCoverage)
        multiplicity = max([len(intercontig.links[0]), len(intercontig.links[1]), multiplicity])
        
        #then, how many times it should appear in each path -> compute the score ?
        totalInteraction = np.sum(interaction_with_path)
        for p in intermediaryContigs[intercontig].keys() :
            estimate = round(multiplicity*interaction_with_path[p]/totalInteraction) #the number of times we expect to find this intermediary in this path
            score += interaction_with_path[p] * 1/(1+abs( intermediaryContigs[intercontig][p] - estimate))
            
    return score
        
#input : untangled paths and the graph
#output : the updated graph with the paths integrated
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
            

        
    
    
    
    
    
    
    
