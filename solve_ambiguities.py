#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

@author: zaltabar

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import matplotlib.pyplot as plt
import numpy as np
import basic_functions as bf
import random
import scipy.integrate as integrate

from shutil import copyfile
from copy import deepcopy 
from transform_gfa import check_links


# this function measures the intensity of interactions between one supercontig and several candidate, including without taking account of the common parts of the supercontigs
# It also weighs the interaction with the length of a supercontig, so that a very long candidate spercontig is not seen as having a lot of connexion just because it is long
def intensity_of_interactions(supercontig, endOfContig, links, candidatesSuperContigs, listOfSuperContigs, interactionMatrix, lengthOfContigs, dist_law, copiesnumber, supercontigsaretouching = True):
    commonContigs = []
    for contig in candidatesSuperContigs[0]:
        common = True
        for sg in candidatesSuperContigs[1:]:
            common = common and (contig in sg)
        if common:
            commonContigs += [contig]
    # now we have a list of common contigs
    for sc in candidatesSuperContigs :
        if len(commonContigs)==len(sc): #that means that the algorithm should wait before taking a decision
            return [-1],[-1]
            
    bestSignature = np.min([copiesnumber[x] for x in supercontig])
    # we return for each supercontig its absolute score and its relative score (wihtout the common parts)
    absoluteScore = [0 for i in range(len(candidatesSuperContigs))]
    relativeScore = [0 for i in range(len(candidatesSuperContigs))]

    total_area = integrate.quad(dist_law, 0, 1000000)[0]
    
    for sg in range(len(candidatesSuperContigs)):
        
        #lengthForNow is a variable keeping track of how far the examined contig of the listOfSuperContigs is from supercontig
        partial_area = 0
        orientation = -1
        for i in links[endOfContig] :
            if listOfSuperContigs[int(i/2)] == candidatesSuperContigs[sg]:
                orientation = i%2 #if supercontig is directly linked to the candidates, then this variable tells us by which end
        lengthForNow = orientation * np.sum([lengthOfContigs[i] for i in candidatesSuperContigs[sg]])
        
        for c in candidatesSuperContigs[sg]:
            
            newLengthForNow = lengthForNow + lengthOfContigs[c] * (0.5-orientation)*2
            if c not in commonContigs :
                #computing the partial area : that way, small supercontigs are not penalized when compared to much longer ones
                partial_area += np.abs(integrate.quad(dist_law, newLengthForNow, lengthForNow)[0])
                #print('Computing partial area for sg ', int(endOfContig/2), ' between ',  newLengthForNow, ' and ', lengthForNow)

            for contig in supercontig:
                if c not in commonContigs and copiesnumber[contig] == bestSignature :
                    absoluteScore[sg] += interactionMatrix[c][contig]
                    relativeScore[sg] += interactionMatrix[c][contig]
                else:
                    absoluteScore[sg] += interactionMatrix[c][contig]
                
            lengthForNow = newLengthForNow
        
        if supercontigsaretouching :
            #print('partial area of contig ', candidatesSuperContigs[sg], ' : ', partial_area)
            relativeScore[sg] *= total_area/partial_area

    return absoluteScore, relativeScore
    
#here we look specifically at one contig and its immediate surroundings (can return -1 if fails in short loop)
def solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs, copiesnumber):
    #we change copiesnumber
    for i in listOfSuperContigs[int(endOfSuperContig/2)]:
        copiesnumber[i] += len(links[endOfSuperContig])-1
    
    #we add all the new supercontigs
    for i in links[endOfSuperContig]:

        otherEndNeighbor = i + 1 - 2 * (i % 2)
        otherEnd = endOfSuperContig + 1 - 2 * (endOfSuperContig % 2)

        if endOfSuperContig in links[otherEndNeighbor] or int(endOfSuperContig/2) == int(i/2) :
            print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(i) +' . Please check that there is indeed a loop there.')
            return -1, -1,-1
        
        if endOfSuperContig%2 == 1 :
            if i%2 == 0 :
                listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig / 2)]+ listOfSuperContigs[int(i / 2)]]
            else :
                listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig / 2)]+ listOfSuperContigs[int(i / 2)][::-1]]
            
            links += [links[otherEnd]] 
            for j in links[otherEnd]:
                links[j] += [len(links) - 1]
            
            links += [links[otherEndNeighbor]]        
            for j in links[otherEndNeighbor]:
                links[j] += [len(links) - 1]
        
        else : 
            if i%2 == 0 :
                listOfSuperContigs += [listOfSuperContigs[int(i / 2)][::-1]+ listOfSuperContigs[int(endOfSuperContig / 2)]]
            else :
                listOfSuperContigs += [listOfSuperContigs[int(i / 2)]+ listOfSuperContigs[int(endOfSuperContig / 2)]]
            
            links += [links[otherEndNeighbor]]             
            for j in links[otherEndNeighbor]:
                links[j] += [len(links) - 1]
            
            links += [links[otherEnd]] 
            for j in links[otherEnd]:
                links[j] += [len(links) - 1]
        
        if endOfSuperContig+1-(endOfSuperContig%2)*2 in links[otherEnd] : #this loop is too difficult for us
            print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(i) +' . Please check that there is indeed a loop there.')    
            return -1,-1,-1
     
    #now we delete the merged supercontigs
    #we start by deleting the links that linked the merged supercontigs to the outside
       
    deletedContigs = [int(endOfSuperContig/2)]
    otherEnd = endOfSuperContig + 1- 2*(endOfSuperContig%2)
    
    #print(endOfSuperContig, listOfSuperContigs[-1], listOfSuperContigs[-2], listOfSuperContigs[int(endOfSuperContig/2)], listOfSuperContigs[int(links[endOfSuperContig][0]/2)], links[otherEnd])
    
    for i in links[otherEnd]:
        try :
            links[i].remove(otherEnd)
        except ValueError : #that means we're in a small loop which whe can't solve
            print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(i) +' . Please check that there is indeed a loop there.')
            return -1, -1,-1
    
    for merged in links[endOfSuperContig]:       
        if len(links[merged]) == 1 : #then the original copy is fully integrated in the supercontig
            deletedContigs += [int(merged/2)]
            otherEnd = merged + 1- 2*(merged%2)
            for i in links[otherEnd]:
                try :
                    links[i].remove(otherEnd)
                except ValueError : #that means we're in a small loop which whe can't solve
                    print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(merged) +' . Please check that there is indeed a loop there.')
                    return -1, -1,-1
        else : #then the original contig still exists by itself, we just delete the link
            try :
                links[merged].remove(endOfSuperContig)
            except ValueError : #that means we're in a small loop which whe can't solve
                    print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(merged) +' . Please check that there is indeed a loop there.')
                    return -1,-1,-1
    
    # then we replace the merged supercontigs and all their links by empty lists (we do not delete them to keep the indexes right)

    for i in deletedContigs:
        listOfSuperContigs[i] = [-1]
        links[2 * i] = [-1]
        links[2 * i + 1] = [-1]
    
    return links, listOfSuperContigs, copiesnumber

# similar to the function above, but simpler : put in one supercontig two smaller supercontig linked by a link unambinguous at both ends
def merge_simply_two_adjacent_contig(endOfSuperContig, links, listOfSuperContigs):

    #print('We are simply putting together ', listOfSuperContigs[int(endOfSuperContig/2)], ' and ', listOfSuperContigs[int(links[endOfSuperContig][0]/2)])
    otherEnd = endOfSuperContig+1-2*(endOfSuperContig%2)
    if links[endOfSuperContig][0] == otherEnd : #then we do not merge a contig with itself
        return -1, -1
    
    otherEndNeighbor = links[endOfSuperContig][0]+1-2*(links[endOfSuperContig][0]%2)
    # we add the new supercontigs
    deleted = links[endOfSuperContig][0]

    if endOfSuperContig%2 == 1 :
             if deleted%2 == 0 :
                 listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig / 2)]+ listOfSuperContigs[int(deleted / 2)]]
             else :
                 listOfSuperContigs += [listOfSuperContigs[int(endOfSuperContig / 2)]+ listOfSuperContigs[int(deleted / 2)][::-1]]
             
             links += [links[otherEnd]] 
             for j in links[otherEnd]:
                 links[j] += [len(links) - 1]
             
             links += [links[otherEndNeighbor]]        
             for j in links[otherEndNeighbor]:
                 links[j] += [len(links) - 1]
         
    else : 
        if deleted%2 == 0 :
            listOfSuperContigs += [listOfSuperContigs[int(deleted / 2)][::-1]+ listOfSuperContigs[int(endOfSuperContig / 2)]]
        else :
            listOfSuperContigs += [listOfSuperContigs[int(deleted / 2)]+ listOfSuperContigs[int(endOfSuperContig / 2)]]
        
        links += [links[otherEndNeighbor]]             
        for j in links[otherEndNeighbor]:
            links[j] += [len(links) - 1]
        
        links += [links[otherEnd]] 
        for j in links[otherEnd]:
            links[j] += [len(links) - 1]

    if endOfSuperContig+1-(endOfSuperContig%2)*2 in links[otherEnd] : #this loop is too difficult for us
        print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(deleted) +' . Please check that there is indeed a loop there.')    
        return -1,-1
    
    if links[otherEnd] == [-1]:
        print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(deleted) +' . Please check that there is indeed a loop there.')    
        return -1,-1

    # now we delete the merged supercontig
    deletedContigs = [int(endOfSuperContig / 2), int(deleted / 2)]

    for i in links[otherEnd]:
        try :
            links[i].remove(otherEnd)
        except ValueError : #that means we're in a small loop with newly created contigs
            print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(i) +' . Please check that there is indeed a loop there.')
            return -1, -1

    for i in links[otherEndNeighbor]:
        links[i].remove(otherEndNeighbor)

    # then we replace the merged supercontigs and all their links by empty lists (we do not delete them for now to keep the indexes right)
    for i in deletedContigs:
        listOfSuperContigs[i] = [-1]
        links[2 * i] = [-1]
        links[2 * i + 1] = [-1]

    return links, listOfSuperContigs

def get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix, lengthOfContigs, dist_law, copiesnumber, thresholdRejected, thresholdAccepted):

    freezed = [] #that's the places where we won't make choices
    
    #loop through all end of contigs, to inspect the robustness of all the links.
    for endOfContig in range(len(links)): 
        if len(links[endOfContig]) > 1: # When there is a 'choice', measurement of the HiC intensity of links
            
            mustBeDeleted = []
            #links are pairwise compared
            
            #first, comparison pairwise the links, those that should be deleted are stored in mustBeDeleted
            for neighbor1 in range(len(links[endOfContig])-1) :
                for neighbor2 in range(neighbor1+1, len(links[endOfContig])):
                    absoluteLinksStrength, linksStrength = intensity_of_interactions(listOfSuperContigs[int(endOfContig / 2)], endOfContig, links, 
                                                                                     [listOfSuperContigs[int(links[endOfContig][neighbor1] / 2)], listOfSuperContigs[int(links[endOfContig][neighbor2] / 2)]],\
                                                                                         listOfSuperContigs,interactionMatrix, lengthOfContigs, dist_law, copiesnumber, True)
                    if absoluteLinksStrength == [-1]: 
                        freezed += [endOfContig]
                        freezed += [eoc for eoc in links[endOfContig]]
                    else :
                        print('I am ', listOfSuperContigs[int(endOfContig/2)],' and here is the choice I have to make : ', [listOfSuperContigs[int(links[endOfContig][neighbor1] / 2)], listOfSuperContigs[int(links[endOfContig][neighbor2] / 2)]], linksStrength)

                        if linksStrength[0]>linksStrength[1] :
                            #     file = open('ratio.txt','a')
                            #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                            if linksStrength[1] < linksStrength[0] * thresholdRejected:  # then that the link does not exist
                                mustBeDeleted.append(neighbor2)
                            elif linksStrength[1] < linksStrength[0] * thresholdAccepted : # then it's not clear, the link is freezed
                                freezed += [endOfContig]
                                
                        else :
                            #     file = open('ratio.txt','a')
                            #     file.write(str(linksStrength[i]/maxStrength)+'\n')
                            if linksStrength[0] < linksStrength[1] * thresholdRejected:  # then that the link does not exist
                                mustBeDeleted.append(neighbor1)
                            elif linksStrength[0] < linksStrength[1] * thresholdAccepted : # then it's not clear, the link is freezed
                                freezed += [endOfContig]
                            
            #second, the links that should be deleted are deleted                
            for i in range(len(links[endOfContig]) - 1, -1, -1):
                if i in mustBeDeleted:
                    print('------Deleting link going from ', listOfSuperContigs[int(endOfContig/2)], ' to ', \
                          listOfSuperContigs[int(links[endOfContig][i]/2)])
                    links[links[endOfContig][i]].remove(endOfContig)
                    del links[endOfContig][i]
            
    return links, freezed

# this function is needed because solve_ambiguity_at_this_end_of_contig creates a lot of useless lists ([-1]) in listOfSuperContigs and in lists, which this fucntion clean
def clean_listOfSuperContigs(links, listOfSuperContigs):

    if links[0] == [-1]:
        list_of_number_of_bad_lists = [1]
    else:
        list_of_number_of_bad_lists = [0]

    for i in links[1:]:
        if i == [-1]:
            list_of_number_of_bad_lists += [list_of_number_of_bad_lists[-1] + 1]
        else:
            list_of_number_of_bad_lists += [list_of_number_of_bad_lists[-1]]

    # now we get rid of all bad lines
    #print('length of links for now : ', len(links))
    #print('length of l')
    links = [i for i in links if i != [-1]]
    listOfSuperContigs = [i for i in listOfSuperContigs if i != [-1]]
    
    for i in range(len(links)):
        for j in range(len(links[i])):
            #print(len(list_of_number_of_bad_lists), links[i][j])
            links[i][j] = links[i][j] - \
                list_of_number_of_bad_lists[links[i][j]]

    return links, listOfSuperContigs

def merge_contigs(links, listOfSuperContigs, copiesnumber, freezed):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends)
    locked = [False for i in listOfSuperContigs]
    
    for endOfSuperContig in range(len(locked)*2):
        #print('end of supercontig that is looked at : ', endOfSuperContig)
        if len(links[endOfSuperContig]) > 1 and endOfSuperContig not in freezed :
            startMerging = not locked[int(endOfSuperContig/2)]
            for i in links[endOfSuperContig] :
                if i < len(locked)*2 :
                    startMerging = startMerging and (not(locked[int(i/2)]))
                else : #meaning we're looking at a freshly created supercontig 
                    startMerging = False
            
            if startMerging : #if nothing is locked for now
                li, lsc, cn = solve_ambiguity_around_this_end_of_contig(endOfSuperContig, deepcopy(links), deepcopy(listOfSuperContigs), copiesnumber.copy())
                if li != -1 : #the function returns -1 when it fell on an unsolvable loop

                    locked[int(endOfSuperContig/2)] = True
                    for i in links[endOfSuperContig] :
                        locked[int(i/2)] = True
                    
                    links, listOfSuperContigs, copiesnumber = deepcopy(li), deepcopy(lsc), deepcopy(cn) #doing all these deepcopies takes time, we might want to reconsider to gain time
                    links = [[links[i][j] for j in range(len(links[i]))] for i in range(len(links))] #we might want to see where that comes from if we want speed
     

    #now we are just going to merge two contigs that are next to each other
    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    
    links, listOfSuperContigs = merge_adjacent_contigs(links, listOfSuperContigs)
    
    return links, listOfSuperContigs, copiesnumber

def merge_adjacent_contigs(links, listOfSuperContigs):
    numberOfInitialLinks = len(links)
    for endOfSuperContig in range(numberOfInitialLinks): 
        if len(links[endOfSuperContig]) == 1 and len(links[links[endOfSuperContig][0]])==1 : #then we just merge
            if links[endOfSuperContig] != [-1] :
                li, lsc = merge_simply_two_adjacent_contig(endOfSuperContig, deepcopy(links), deepcopy(listOfSuperContigs))  
                if li != -1 :
                    links, listOfSuperContigs = deepcopy(li), deepcopy(lsc)
                    links = [[links[i][j] for j in range(len(links[i]))] for i in range(len(links))]

    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    return links, listOfSuperContigs
    
def solve_ambiguities(links, names, interactionMatrix, lengthOfContigs, dist_law, stringenceReject, stringenceAccept, steps):

    copiesnumber = [1 for i in names]
    listOfSuperContigs = [[x] for x in range(len(names))] #The contigs are numbered in listOfSuperContigs, correspondance can be made in the list 'names'
    
    links, listOfSuperContigs = merge_adjacent_contigs(links, listOfSuperContigs)
    
    for i in range(steps):
        
        print(str(i / steps * 100) + "% of solving ambiguities done\n")

        links, freezed = get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix, lengthOfContigs, dist_law, copiesnumber, stringenceReject, stringenceAccept)

        links, listOfSuperContigs, copiesnumber = merge_contigs(links, listOfSuperContigs, copiesnumber, freezed)
        bf.export_to_GFA(links, listOfSuperContigs, copiesnumber, names, exportFile = 'tests/fake'+str(i)+'.gfa')
        
    return links, listOfSuperContigs, copiesnumber

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


