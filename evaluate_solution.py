#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:37:22 2020

File dedicated to functions evaluating the quality of the final GFA
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from copy import deepcopy

import input_output as io

from input_output import load_gfa
from solve_ambiguities import solve_ambiguities

#function scoring the output of solve_ambiguities (evaluating how well HiC contacts correlate with GFA distance)
def score_output(listOfSuperContigs, links, lengthOfContigs, interactionMatrix, infinite_distance = 500000):
    
    #fist, determine the shortest distance in the graph between two contigs
    matrixOfShortestPaths = find_matrixOfShortestPaths(listOfSuperContigs, links, lengthOfContigs, infinite_distance)
    
    #Now computing the score : the higher the score, the less correlation there is between HiC contacts and distance in GFA : the goal is to build a GFA with a low score
    score = 0
    
    n = len(lengthOfContigs)
    interactionMatrixMean = np.sum([np.sum([interactionMatrix[i][j] for j in range(i+1, n)]) for i in range(n-1)])*2/n/(n-1)
    
    #the first element of score is a sum of distance*HiC contacts, encouraging contigs with big HiC contacts to be close
    for i in range(len(lengthOfContigs)-1):
            for j in range(i+1, len(lengthOfContigs)) :
                score += (interactionMatrix[i][j]-interactionMatrixMean) * matrixOfShortestPaths[i][j]/infinite_distance/interactionMatrixMean #-interactionMatrixMean so that two contigs that do not have HiC contacts are encouraged to be far away
    
    score /= n*(n-1)/2 # that ensures that score cannot be smaller than -1
    #the second element of score is a parcimonious element : duplication of contigs is penalized, as to not duplicate freely the contigs to minimize the first term
    # for c in range(len(copiesnumber)) :
    #     score += (copiesnumber[c]-1)*infinite_distance*np.max(interactionMatrix[i])*0.1

    return score
    
def find_matrixOfShortestPaths(listOfSuperContigs, links, lengthOfContigs, infinite_distance = 500000): #infinite_distance is a int parameter : if two contigs don't touch, they will be infinite_distance apart
    
    #Initializing by building a "matrix" in which all elements of listOfSuperContigs have a distance with all elements of listOfSuperContigs
    matrixOfShortestPathWithContigsReplication = []
    for i in range(len(listOfSuperContigs)) :
        matrixOfShortestPathWithContigsReplication.append([])
        for j in range(len(listOfSuperContigs[i])) :
            matrixOfShortestPathWithContigsReplication[-1].append([[infinite_distance for k in l] for l in listOfSuperContigs])
            matrixOfShortestPathWithContigsReplication[-1][-1][i][j] = -lengthOfContigs[listOfSuperContigs[i][j]]
            
     
    #Now running Dijkstra algorithm : its original aim is to find the shortest path between two points. Here, by memorizing distances, it is generalized at finding the shortest distance between any two given contigs
    Continue = True
    while(Continue) :
        Continue = False
        
        for i in range(len(listOfSuperContigs)):
            for j in range(len(listOfSuperContigs[i])):
                
                if len(listOfSuperContigs[i]) == 1 :
                    
                    for touchingSuperContig in links[i*2]:

                        for k in range(len(listOfSuperContigs)):
                            for l in range(len(listOfSuperContigs[k])):
                                if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                    matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]] :
                                        Continue = True
                                        matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]]

                    for touchingSuperContig in links[i*2+1]:

                        for k in range(len(listOfSuperContigs)):
                            for l in range(len(listOfSuperContigs[k])):
                                if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                    matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]] :
                                        Continue = True
                                        matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]]

                
                elif j == 0 : #looking at the beginning of a supercontig
          
                    for touchingSuperContig in links[i*2]:

                        for k in range(len(listOfSuperContigs)):
                            for l in range(len(listOfSuperContigs[k])):
                                if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                    matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]] :
                                        Continue = True
                                        matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]]
    
                    for k in range(len(listOfSuperContigs)):
                        for l in range(len(listOfSuperContigs[k])):
                            if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                matrixOfShortestPathWithContigsReplication[i][j+1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]] :
                                    Continue = True
                                    matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[i][j+1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]]    
                
                elif j == len(listOfSuperContigs[i])-1 : #looking at the end of a supercontig
                
                    for touchingSuperContig in links[i*2+1]:

                        for k in range(len(listOfSuperContigs)):
                            for l in range(len(listOfSuperContigs[k])):
                                if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                    matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]] :
                                        Continue = True
                                        matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[int(touchingSuperContig/2)][-(touchingSuperContig%2)][k][l] + lengthOfContigs[listOfSuperContigs[int(touchingSuperContig/2)][-(touchingSuperContig%2)]]
    
                    for k in range(len(listOfSuperContigs)):
                        for l in range(len(listOfSuperContigs[k])):
                            if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                matrixOfShortestPathWithContigsReplication[i][j-1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]] :
                                    Continue = True
                                    matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[i][j-1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]]    
     
                
                else : #looking at the middle of a supercontig

                    for k in range(len(listOfSuperContigs)):
                        for l in range(len(listOfSuperContigs[k])):
                            if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                matrixOfShortestPathWithContigsReplication[i][j-1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]] :
                                    Continue = True
                                    matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[i][j-1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]]
                            if matrixOfShortestPathWithContigsReplication[i][j][k][l] > \
                                matrixOfShortestPathWithContigsReplication[i][j+1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]] :
                                    Continue = True
                                    matrixOfShortestPathWithContigsReplication[i][j][k][l] = matrixOfShortestPathWithContigsReplication[i][j+1][k][l] + lengthOfContigs[listOfSuperContigs[i][j]]
        
    #Now putting together the contigs that were duplicated in listOfSuperContigs 
    matrixOfShortestPath = [[infinite_distance for i in lengthOfContigs] for j in lengthOfContigs]
    
    for i in range(len(listOfSuperContigs)):
        for j in range(len(listOfSuperContigs[i])):
            for k in range(len(listOfSuperContigs)) :
                for l in range(len(listOfSuperContigs[k])) :
                    
                    if matrixOfShortestPath[listOfSuperContigs[i][j]][listOfSuperContigs[k][l]] > matrixOfShortestPathWithContigsReplication[i][j][k][l] :
                        matrixOfShortestPath[listOfSuperContigs[i][j]][listOfSuperContigs[k][l]] = matrixOfShortestPathWithContigsReplication[i][j][k][l]

    #taking out the negative values of the matrix
    for i in matrixOfShortestPath :
        for j in i :
            if j < 0 :
                j = 0
                
    return matrixOfShortestPath

def draw_distance_HiCcontacts_correlation(listOfSuperContigs, links, lengthOfContigs, interactionMatrix, infinite_distance = 500000):
    
    #fist, determine the shortest distance in the graph between two contigs
    matrixOfShortestPaths = find_matrixOfShortestPaths(listOfSuperContigs, links, lengthOfContigs, infinite_distance)
    print('shortest paths found')
    
    shortestPaths = []
    HiCcontacts = []
    
    
    for i in range(len(lengthOfContigs)-1):
            for j in range(i+1, len(lengthOfContigs)) :
                if interactionMatrix[i][j] > 10 :
                    HiCcontacts.append (interactionMatrix[i][j]/lengthOfContigs[i]/lengthOfContigs[j])
                    shortestPaths.append(matrixOfShortestPaths[i][j])
      
    print('coucou')
    plt.scatter(HiCcontacts, shortestPaths, alpha = 0.1)
    plt.xlabel('Intensity of contacts')
    plt.ylabel('Smallest distance between the two contigs')
    plt.show()

def heat_solution(links, listOfSuperContigs, copiesNumber, originalLinks, heat) :
    
    #merge, with a probability depending on the heat, contigs that have been previously duplicated
    for contigToMerge in range(len(copiesNumber)) :
        
        if copiesNumber[contigToMerge] > 1 :
            if random.random() < heat : #then start merging
                
                copiesNumber[contigToMerge] = 1
                
                originalLength = len(listOfSuperContigs)

                listOfSuperContigs.append([contigToMerge])
                links.append([])
                links.append([])
                for supercontig, sc in enumerate(listOfSuperContigs[:originalLength]) :
                    
                    if contigToMerge in sc :
                
                        print(listOfSuperContigs)
                        print(sc, contigToMerge)
                        print(links[supercontig*2+1])
                        print(links)
                        if len(sc) == 1 :
                            
                            links[originalLength*2] += links[supercontig*2]
                            for i in links[supercontig*2]:
                                links[i].append(originalLength*2)
                                links[i].remove(supercontig*2)
                                
                            links[originalLength*2+1] += links[supercontig*2+1]
                            for i in links[supercontig*2+1]:
                                links[i].remove(supercontig*2+1)
                                links[i].append(originalLength*2+1)
                            
                            listOfSuperContigs[supercontig] = [-1]
                            links[supercontig*2] = [-1]
                            links[supercontig*2+1] = [-1]
                            
                        elif contigToMerge == sc[0] :
                            links[originalLength*2] += links[supercontig*2]
                            for i in links[supercontig*2]:
                                links[i].remove(supercontig*2)
                                links[i].append(originalLength*2)
                                
                            listOfSuperContigs.append(sc[1:])

                            links.append([originalLength*2+1])
                            links[originalLength*2+1] += [len(links)-1]

                            links.append(links[supercontig*2+1])
                            for i in links[supercontig*2+1]:
                                links[i].remove(supercontig*2+1)
                                links[i].append(len(links)-1)
                            
                            listOfSuperContigs[supercontig] = [-1]
                            links[supercontig*2] = [-1]
                            links[supercontig*2+1] = [-1]
                            
                        elif contigToMerge == sc[-1] :
                            links[originalLength*2+1] += links[supercontig*2+1]
                            for i in links[supercontig*2+1]:
                                links[i].remove(supercontig*2+1)
                                links[i].append(originalLength*2+1)

                            listOfSuperContigs.append(sc[:len(sc)-1])
                            links.append(links[supercontig*2])
                            for i in links[supercontig*2]:
                                links[i].remove(supercontig*2)
                                links[i].append(len(links)-1)
                            links.append([originalLength*2])
                            links[originalLength*2] += [len(links)-1]
                            
                            listOfSuperContigs[supercontig] = [-1]
                            links[supercontig*2] = [-1]
                            links[supercontig*2+1] = [-1]

                        else :
                            indexOfContig = sc.index(contigToMerge)
                            
                            listOfSuperContigs.append(sc[:indexOfContig])
                            links.append(links[supercontig*2])
                            for i in links[supercontig*2]:
                                links[i].remove(supercontig*2)
                                links[i].append(len(links)-1)
                            links.append([originalLength*2])
                            
                            links[originalLength*2] += [len(links)-1]
                        
                            listOfSuperContigs.append(sc[indexOfContig+1:len(sc)])
                            links.append([originalLength*2+1])
                            links[originalLength*2+1] += [len(links)-1]
                            
                            links.append(links[supercontig*2+1])
                            for i in links[supercontig*2+1]:
                                links[i].remove(supercontig*2+1)
                                links[i].append(originalLength*2+1)
                                                              
                            listOfSuperContigs[supercontig] = [-1]
                            links[supercontig*2] = [-1]
                            links[supercontig*2+1] = [-1]
                            
                        for l in range(len(links)) : #to delete all redundant links we may have created
                            links[l] = list(set(links[l]))

                #print(listOfSuperContigs)
                
                links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
                
    return links, listOfSuperContigs, copiesNumber
    
def simulated_annealing(links, names, interactionMatrix, lengthOfContigs, dist_law, stringenceReject, stringenceAccept, steps) :
    
    originalLinks = deepcopy(links)
    
    #run solve_ambiguities the first time :
    links, listOfSuperContigs, copiesNumber = solve_ambiguities(links, names, interactionMatrix, lengthOfContigs, dist_law, stringenceReject, stringenceAccept, steps)
    best_score = score_output(listOfSuperContigs, links, lengthOfContigs, interactionMatrix)    
    
    for i in range(10) :
        #first, modifiy a bit the solution previously found
        newLinks, newListOfSuperContigs, newCopiesNumber = heat_solution(deepcopy(links), deepcopy(listOfSuperContigs), deepcopy(copiesNumber), originalLinks, heat = (10-i)/20)
        
        #secondly, proceed to solve the new problem
        newLinks, newListOfSuperContigs, newCopiesNumber = solve_ambiguities(newLinks, names, interactionMatrix, lengthOfContigs,dist_law, stringenceReject, stringenceAccept, steps, newListOfSuperContigs, newCopiesNumber)
        
        #thirdly, score the new output to see if it is more satisfying than the last one
        newScore = score_output(newListOfSuperContigs, newLinks, lengthOfContigs, interactionMatrix)  
        
        if newScore < best_score :
            best_score = newScore
            links = newLinks
            listOfSuperContigs = newListOfSuperContigs
            copiesNumber = newCopiesNumber
        
    return links, listOfSuperContigs, copiesNumber
        
# originalLinks, names, lengthOfContigs = load_gfa('results/A_Vaga_finished.gfa')
# #info_contigs = bf.read_info_contig('data/results/info_contigs.txt')
# interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')

# #let's build the new interaction matrix
# newInteractionMatrix = []
# for i in range(len(names)) :
#     newInteractionMatrix.append([])
#     for j in range(len(names)):
#         newInteractionMatrix[-1].append(interactionMatrix[int(int(names[i].split('-')[0])/2)][int(int(names[j].split('-')[0])/2)])
     
# print('Loaded')
# #print(score_output([[i] for i in range(len(names))], originalLinks, lengthOfContigs, newInteractionMatrix, infinite_distance = 1000000))
# draw_distance_HiCcontacts_correlation([[i] for i in range(len(names))], originalLinks, lengthOfContigs, newInteractionMatrix, infinite_distance = 1000000)

# print('Finished')