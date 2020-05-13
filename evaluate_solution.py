#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:37:22 2020

File dedicated to functions evaluating the quality of the final GFA
"""

#function scoring the output of solve_ambiguities (evaluating how well HiC contacts correlate with GFA distance)
def score_output(listOfSuperContigs, links, lengthOfContigs, interactionMatrix):
    
    #fist, determine the shortest distance in the graph between two contigs
    matrixOfShortestPaths = find_matrixOfShortestPaths(listOfSuperContigs, links, lengthOfContigs)
    
    #Now computing the score : the higher the score, the less correlation there is between HiC contacts and distance in GFA : the goal is to build a GFA with a low score
    score = 0
    
    #the first element of score is a sum of distance*HiC contacts, encouraging contigs with big HiC contacts to be close
    for i in range(len(lengthOfContigs)-1):
            for j in range(i+1, len(lengthOfContigs)) :
                score += interactionMatrix[i][j] * matrixOfShortestPaths[i][j]
                
    #the second element of score is a parcimonious element : duplication of contigs is penalized, as to not duplicate needlessly the contigs
    return 0
    
def find_matrixOfShortestPaths(listOfSuperContigs, links, lengthOfContigs, infinite_distance = 500000): #infinite_distance is a int parameter : if two contigs don't touch, they will be infinite_distance apart
    
    #Initializing by building a "matrix" in which all elements of listOfSuperContigs have a distance with all elements of listOfSuperContigs
    matrixOfShortestPathWithContigsReplication = []
    for i in range(len(listOfSuperContigs)) :
        matrixOfShortestPathWithContigsReplication.append([])
        for j in range(len(listOfSuperContigs[i])) :
            matrixOfShortestPathWithContigsReplication[-1].append([[infinite_distance for k in l] for l in listOfSuperContigs])
            matrixOfShortestPathWithContigsReplication[-1][-1][i][j] = -lengthOfContigs[listOfSuperContigs[i][j]]
            
    #print(matrixOfShortestPathWithContigsReplication)
     
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
    return matrixOfShortestPath


print(find_matrixOfShortestPaths([[0,1,2],[0]], [[],[2],[1],[]], [10000 for i in range(3)]))