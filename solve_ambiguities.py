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
from shutil import copyfile

# this function measures the intensity of interactions between one supercontig and several candidate, including without taking account of the common parts of the supercontigs
# the contig in this function are not numbered by their end, i.e. give it [1234] and not [2468,2469]
def intensity_of_interactions(supercontig, listOfSuperContigs, interactionMatrix):
    commonContigs = []
    for contig in listOfSuperContigs[0]:
        common = True
        for sg in listOfSuperContigs[1:]:
            common = common and (contig in sg)
        if common:
            commonContigs += [contig]
    # now we have a list of common contigs

    # we return for each supercontig its absolute score and its relative score (wihtout the common parts)
    absoluteScore = [0 for i in range(len(listOfSuperContigs))]
    relativeScore = [0 for i in range(len(listOfSuperContigs))]

    for sg in range(len(listOfSuperContigs)):
        for c in listOfSuperContigs[sg]:
            for contig in supercontig:
                if c not in commonContigs:
                    absoluteScore[sg] += interactionMatrix[c][contig]
                    relativeScore[sg] += interactionMatrix[c][contig]
                else:
                    absoluteScore[sg] += interactionMatrix[c][contig]

    return absoluteScore, relativeScore
    
#here we look specifically at one contig and its immediate surroundings (can return -1 if fails in short loop)
def solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs, copiesnumber):
            
    #we change copiesnumber
    for i in listOfSuperContigs[int(endOfSuperContig/2)]:
        copiesnumber[i] += len(links[endOfSuperContig])-1
    
    #we add all the new supercontigs
    for i in links[endOfSuperContig]:
        listOfSuperContigs += [
            listOfSuperContigs[int(endOfSuperContig / 2)]
            + listOfSuperContigs[int(i / 2)]
        ]

        otherEnd = endOfSuperContig + 1 - 2 * (endOfSuperContig % 2)
        links += [links[otherEnd]]
        for j in links[otherEnd]:
            links[j] += [len(links) - 1]

        otherEnd = i + 1 - 2 * (i % 2)
        links += [links[otherEnd]]
        
        if endOfSuperContig+1-(endOfSuperContig%2)*2 in links[otherEnd] : #this loop is too difficult for us
            print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(i) +' . Please check that there is indeed a loop there.')    
        return -1,-1,-1
        for j in links[otherEnd] :
            links[j] += [len(links)-1]
            
    #now we delete the merged supercontigs
    #we start by deleting the links that linked the merged supercontigs to the outside
       
    deletedContigs = [int(endOfSuperContig/2)]
    otherEnd = endOfSuperContig + 1- 2*(endOfSuperContig%2)
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

    if links[endOfSuperContig][0] == endOfSuperContig+1-2*(endOfSuperContig%2) : #then we do not merge a contig with itself
        return -1, -1
    
    # we add the new supercontigs
    deleted = links[endOfSuperContig][0]
    listOfSuperContigs += [
        listOfSuperContigs[int(endOfSuperContig / 2)]
        + listOfSuperContigs[int(deleted / 2)]
    ]

    otherEnd = endOfSuperContig + 1 - 2 * (endOfSuperContig % 2)
    links += [links[otherEnd]]

    if endOfSuperContig+1-(endOfSuperContig%2)*2 in links[otherEnd] : #this loop is too difficult for us
        print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(deleted) +' . Please check that there is indeed a loop there.')    
        return -1,-1
    
    for j in links[otherEnd]:
        links[j] += [len(links) - 1]

    otherEnd = deleted + 1 - 2 * (deleted % 2)
    links += [links[otherEnd]]
    
    if links[otherEnd] == [-1]:
        print('We have a difficulty here : '+str(endOfSuperContig)+ ' ' + str(deleted) +' . Please check that there is indeed a loop there.')    
        return -1,-1
    for j in links[otherEnd]:
        links[j] += [len(links) - 1]

    # now we delete the merged supercontig
    deletedContigs = [int(endOfSuperContig / 2), int(deleted / 2)]
    otherEnd = endOfSuperContig + 1 - 2 * (endOfSuperContig % 2)
    for i in links[otherEnd]:
        links[i].remove(otherEnd)

    otherEnd = deleted + 1 - 2 * (deleted % 2)
    for i in links[otherEnd]:
        links[i].remove(otherEnd)

    # then we replace the merged supercontigs and all their links by empty lists (we do not delete them to keep the indexes right)
    for i in deletedContigs:
        listOfSuperContigs[i] = [-1]
        links[2 * i] = [-1]
        links[2 * i + 1] = [-1]

    return links, listOfSuperContigs


def get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix, strengthThreshold=3):

    for endOfContig in range(len(links)):
        if len(links[endOfContig]) > 1:
            absoluteLinksStrength, linksStrength = intensity_of_interactions(
                listOfSuperContigs[int(endOfContig / 2)],
                [listOfSuperContigs[int(x / 2)] for x in links[endOfContig]],
                interactionMatrix,
            )
            maxStrength = np.max(linksStrength)

            for i in range(len(links[endOfContig]) - 1, -1, -1):
                if (
                    linksStrength[i] < maxStrength / strengthThreshold
                ):  # we consider then that the link does not exist
                    links[links[endOfContig][i]].remove(endOfContig)
                    del links[endOfContig][i]
    return links


# this function is needed because solve_ambiguity_at_this_end_of_contig creates a lot of useless lists ([-1]) in listOfSuperContigs and in lists
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
    links = [i for i in links if i != [-1]]
    listOfSuperContigs = [i for i in listOfSuperContigs if i != [-1]]

    for i in range(len(links)):
        for j in range(len(links[i])):
            links[i][j] = links[i][j] - list_of_number_of_bad_lists[links[i][j]]

    return links, listOfSuperContigs

def merge_contigs(links, listOfSuperContigs, copiesnumber):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends)
    locked = [False for i in listOfSuperContigs]
    
    for endOfSuperContig in range(len(locked)*2):
        #print('end of supercontig that is looked at : ', endOfSuperContig)
        if len(links[endOfSuperContig]) > 1 :
            startMerging = not locked[int(endOfSuperContig/2)]
            for i in links[endOfSuperContig] :
                if i < len(locked)*2 :
                    startMerging = startMerging and (not(locked[int(i/2)]))
                else :
                    startMerging = False
            
            if startMerging : #if nothing is locked for now
                li, lsc, cn = solve_ambiguity_around_this_end_of_contig(endOfSuperContig, links, listOfSuperContigs, copiesnumber)
                if li != -1 : #the function returns -1 when it fell on an unsolvable loop
                    links, listOfSuperContigs, copiesnumber = li, lsc, cn
                    links = [[links[i][j] for j in range(len(links[i]))] for i in range(len(links))] #we might want to see where that comes from if we want speed
                    locked[int(endOfSuperContig/2)] = True
                    for i in links[endOfSuperContig] :
                        locked[int(i/2)] = True
            
    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    locked = [False for i in listOfSuperContigs]

    # for endOfSuperContig in range(len(locked)*2):
    #     if len(links[endOfSuperContig]) == 1 and len(links[links[endOfSuperContig][0]])==1 : #then we just merge
    #         if links[endOfSuperContig] != [-1] :
    #             li, lsc = merge_simply_two_adjacent_contig(endOfSuperContig, links, listOfSuperContigs)  
    #             if li != -1 :
    #                 links, listOfSuperContigs = li, lsc
    #                 links = [[links[i][j] for j in range(len(links[i]))] for i in range(len(links))]

    links, listOfSuperContigs = clean_listOfSuperContigs(links, listOfSuperContigs)
    return links, listOfSuperContigs, copiesnumber

def export_to_GFA(links, listOfSuperContigs, originalLinks, copiesnumber, fastaFile):

  #  copyfile("results/sequencesWithoutLinks.gfa", "results/newAssembly.gfa")
    f = open("results/newAssembly.gfa", "w")
    addresses = [-1 for i in range(len(listOfSuperContigs)*2)]
    copiesUsed = [0 for i in copiesnumber]
    
    for sc,supercontig in enumerate(listOfSuperContigs) :
        if sc%30 == 0 :
            print(int(sc/len(listOfSuperContigs)*1000)/10, '% of exporting done')
        for c, contig in enumerate(supercontig) :
            f.write('S\t'+str(contig*2)+'-'+str(copiesUsed[contig])+'\t'+bf.get_contig(fastaFile, supercontig[0], supercontig[0]*2-1)+'\n')
            #print(supercontig)
            if c > 0:
                print('coucou')
                f.write('L\t'+str(supercontig[c-1]*2)+'-'+str(copiesUsed[supercontig[c-1]]-1)+\
                        '\t+\t' + str(contig*2)+'-'+\
                            str(copiesUsed[contig])+'\t+\t*\n')
            
            copiesUsed[contig] += 1


def solve_ambiguities(links, listOfContigs, interactionMatrix, copiesnumber):  # look at ambilguities one after the other
    listOfSuperContigs = [[x] for x in listOfContigs]
    steps = 1
    for i in range(steps):
        if i % 100 == 0:
            print(str(i / steps * 100) + "% of solving ambiguities done")

        links = get_rid_of_bad_links(links, listOfSuperContigs, interactionMatrix)
        links, listOfSuperContigs, copiesnumber = merge_contigs(links, listOfSuperContigs, copiesnumber)

    return links, listOfSuperContigs, copiesnumber

links = bf.import_links('listsPython/links.csv')
# print(links[1465])
# infContigs = bf.read_info_contig('data/results/info_contigs.txt')
interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')
for i in range(len(interactionMatrix)):
     interactionMatrix[i][i] = 0
print('Loaded')

# links, listOfSuperContigs = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix)
# print(listOfSuperContigs)
# print('Now the end of links ')
# print(links[2600:])
# print(intensity_of_interactions([874*2, 874*2+1], [[584*2,584*2+1,1120*2,1120*2+1], [584*2,584*2+1,78*2,78*2+1]], interactionMatrix))
# print(intensity_of_interactions([584*2,584*2+1,78*2,78*2+1], [[802*2,802*2+1,874*2, 874*2+1], [802*2,802*2+1,743*2,743*2+1]], interactionMatrix))

# print(
#     solve_ambiguities(
#         [[], [4], [], [4], [1, 3], [6, 8], [5], [], [5], []],
#         [0, 1, 2, 3, 4],
#         [
#             [1, 1, 1, 1, 0],
#             [1, 1, 1, 0, 1],
#             [1, 1, 1, 1, 1],
#             [1, 0, 1, 1, 1],
#             [0, 1, 1, 1, 1],
#         ],
#     )
# )
        
newlinks, listOfSuperContigs, copiesnumber = solve_ambiguities(links, [x for x in range(1312)], interactionMatrix, [1 for i in range(1312)])
print('Now exporting')
print(listOfSuperContigs)
#print(copiesnumber)
export_to_GFA(newlinks, listOfSuperContigs, links, copiesnumber, 'data/Assembly.fasta') 
                        
print('Finished')


