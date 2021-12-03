#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 18:06:07 2021

@author: rfaure

A file summarizing the algorithmic part of untangling the graph with HiC
"""

from determine_multiplicity import determine_multiplicity
from scipy import sparse
from scipy.sparse import isspmatrix_csr
import matplotlib.pyplot as plt

#input : segments and the interactionMatrix of Hi-C contacts. Optionnaly, a list of haploid contigs obtained from the long reads algorithm.
def solve_with_HiC(segments, interactionMatrix, names):
    
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
    contacts = sparse.lil_matrix((len(haploidContigs)*2, len(haploidContigs)*2))
    sure = False
    #while not sure :
    haploidContigs, haploidContigsNames, sure =  match_haploidContigs(segments, names, normalInteractions, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames)


#input: a graph, in the form of a list of segments and a list of haploid conigs
#output : a list of all independant knots of the graph. A new list of haploidContigs where each haploid contig actually has contacts with neighboring haploid contigs
def determine_list_of_knots(segments, multiplicities, haploidContigs, haploidContigsNames, interactionMatrix, names) :
    
    #first compute the list of neighbors for all haploidContigs, checking if each contig is informative or not
    listOfNeighbors = [[] for i in range(len(haploidContigs)*2)]
    notInformative = []
    for end in range(len(listOfNeighbors)) :
        segmentsAlreadyTraversed = set()
        listOfNeighbors[end] = find_neighbors(haploidContigs[end//2], end%2, segmentsAlreadyTraversed, haploidContigsNames)
        
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
            
            # if 'edge_358' in haploidContigs[end//2].names :
            #     print("That's edge 358 : ", interactions, " ", [haploidContigs[i//2].names for i in listOfNeighbors[end]])
            #     while True :
            #         r=0
            if all([i==0 for i in interactions]) :
                notInformative += [end//2]
            
    #remove all uninformative haploidContigs from haploidContigs
    print("Here are all the uninformative contigs I should remove : ", [haploidContigs[i].names for i in notInformative])
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

    #Now move on to the untangling
    
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
def match_haploidContigs(segments, names, interactionMatrix, list_of_neighbors, list_of_knots, contacts, knotOfContig, haploidContigs, haploidContigsNames) :
    
    sure_haploids = True
    not_actually_haploid = [] #a list of contigs that have no contacts with neighbors among the haploidContigs (those contigs are useless)
    
    for knot in list_of_knots :
        
        if len(knot) > 1 :
            preferred_contact = {}
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
                    #interactions += [interactionMatrix[index][ names[haploidContigs[neighbor//2].full_name()]]]
                #interactions = [ interactionMatrix[index][ haploidContigsNames[haploidContigs[i//2].full_name()] ] for i in knot ]
                print("Looking at interaction from contig ", haploidContigs[end//2].full_name(), " and here are its interactions: ", interactions)
                
                m = max(interactions)
                
                if m > 0 :
                    bestIdx = interactions.index(max(interactions))
                    
                    if list_of_neighbors[end][bestIdx] not in preferred_contact :
                        preferred_contact[list_of_neighbors[end][bestIdx]] = end  
                    else : #it means two contigs point to the same contig, which is a problem
                        #knot_solved = False
                        #print("Contig ", haploidContigs[list_of_neighbors[end][bestIdx]//2].names, " does not look so haploid")
                        not_actually_haploid += [list_of_neighbors[end][bestIdx]//2]
                        sure_haploids = False
                    
                    contacts[list_of_neighbors[end][bestIdx], end] = 1
                    contacts[end, list_of_neighbors[end][bestIdx]] = 1
                    
                else :
                    print("I can't do anything about contig ", haploidContigs[end//2].names)
                    knot_solved = False

                    
            if knot_solved :
                print("We solved a knot:")
                for e, end in enumerate(knot) :
                    for otherEnd in range(contacts.shape[0]) :
                        if contacts[otherEnd, end] == 1 :
                            print(haploidContigs[end//2].names, " -> ", haploidContigs[otherEnd//2].names)
            print()
         
    #now get rid of non-haploid contigs
    
    #first update listOfKnots, list_of_neighbors and knotOfContig with the new information on haploidContigs
    
    # for not_hap in not_actually_haploid :
    #     #update knotOfContig
    #     if 

    #then update the haploidContigs list
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
            
    return reliable_haploid_contigs, reliable_haploid_contigsNames, sure_haploids #return the updated list of haploid contigs

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




