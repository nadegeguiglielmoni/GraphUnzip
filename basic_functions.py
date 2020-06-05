#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""
import pandas as pd
import numpy as np
# import sparse module from SciPy package 
from scipy import sparse #dok_matrix will be used for interactionMatrix since they allow for fast access of elements

# NOT USED
# Read sparse matrix
# Input :
#   file : sparse matrix with 3 fields: frag1, frag2, contacts
# Output :
#   content : list of parsed contacts
def read_abs_fragments_contact_weighted(file, header=True):

    content = []

    with open(file) as f:
        for line in f:
            content.append([line])

    # parsing lines
    content = [x.strip("\n") for x in content]
    content = [x.split("\t") for x in content]
    content = [[int(x[0]), int(x[1]), int(x[2])] for x in content]

    # remove header
    if header:
        content = content[1:]

    return content


# Read fragments list file
# /!\ for many functions to work, fragmentList has to be sorted (first contig1, then contig2...)
# Input :
#   file : fragments_list.txt
# Output :
#   content : list with contig_id, fragment start, fragment end, fragment length, where the list
#             index is the fragment id
def read_fragment_list(file, header=True):

    with open(file) as f:
        content = f.readlines()

    if header:
        del content[0]

    # parsing
    # 1: contig_id, 2: fragment_start, 3: fragment_end, 4: fragment_length
    content = [x.strip("\n").split("\t") for x in content]
    content = [
        [x[1].strip("sequence"), int(x[2]), int(x[3]), int(x[4])] for x in content
    ]
    # the removal of "sequence" is too specific to a certain contig name format
    # should be removed/improved

    return content


def read_info_contig(file):

    with open(file) as f:
        content = f.readlines()

    # parsing
    content = [x.strip("\n") for x in content[1:]]
    content = [x.strip("sequence") for x in content]
    content = [x.split("\t") for x in content]
    content = [[int(x[0]), int(x[1]), int(x[2]), int(x[3])] for x in content]

    return content


def interactionMatrix(hiccontactsfile, fragmentList, header=True):  # the header refers to the hiccontactsfile

    # create a full interaction matrix of contig vs contig
    # 1 -> [1...N] N contigs
    # ...
    # N -> [1...N]
    print('Starting to build an interaction matrix')
    interactionMatrix = sparse.dok_matrix((int(fragmentList[-1][0]) + 1, int(fragmentList[-1][0]) + 1))
    # interactionMatrix = [
    #     [0 for i in range(int(fragmentList[-1][0]) + 1)]
    #     for j in range(int(fragmentList[-1][0]) + 1)]
    
    n = 0
    with open(hiccontactsfile) as f:
        
        for line in f:
    
            if not header :
                line = line.strip("\n").split("\t")
        
                # frag1, frag2, contacts
                contact = [int(line[0]), int(line[1]), int(line[2])]
        
                # search for contig name corresponding to fragment id
                contig1 = int(fragmentList[contact[0]][0])
                contig2 = int(fragmentList[contact[1]][0])
        
                # add contacts to interaction matrix
                interactionMatrix[contig1, contig2] += contact[2]
                interactionMatrix[contig2, contig1] += contact[2]
                
                n += 1
                if n%10000 == 0 :
                    print(n/10000, '*10000 hic contacts inputed in the matrix')
            else :
                header = False

    return interactionMatrix


def export_to_csv(l, file):
    df = pd.DataFrame(l)
    df.to_csv(file)


def import_from_csv(file):

    df = pd.read_csv(file)
    l = df.values.tolist()

    newl = [[x for x in i if not pd.isnull(x)] for i in l]
    return [x[1:] for x in newl]  # we discard the header line


def import_links(file):
    links = import_from_csv(file)
    return [[int(i) for i in j] for j in links]


def get_contig(fastaFile, contig, firstline=0):

    with open(fastaFile) as f:

        lookAtNextLine = False
        linenumber = 0
        for line in f:
            if linenumber >= firstline:
                if lookAtNextLine:
                    return line
                target = ">" + str(contig)
                if target in line:
                    lookAtNextLine = True

            linenumber += 1
    return "In get_contig : the contig you are seeking is not in the fasta file"

def get_contig_GFA(gfaFile, contig):
    
    with open(gfaFile) as f:

        for line in f:
            
            sline = line.split('\t')
            if len(sline) >= 3 :
                if sline[0] == 'S' and (contig in sline[1]) :
                    return sline[2]

    return "In get_contig : the contig you are seeking is not in the fasta file"

def export_to_GFA(
    links,
    listOfSuperContigs,
    copiesnumber,
    originalLinks,
    originalLinksCIGAR = None,
    names=None,
    gfaFile="",
    exportFile="results/newAssembly.gfa",
):

    if names == None:
        names = [str(2 * i) for i in range(int(len(links) / 2))]
    #  copyfile("results/sequencesWithoutLinks.gfa", "results/newAssembly.gfa")
    f = open(exportFile, "w")
    addresses = ["" for i in range(len(listOfSuperContigs) * 2)]
    copiesUsed = [0 for i in copiesnumber]

    #fist write the sequences and the links within the supercontigs
    for sc, supercontig in enumerate(listOfSuperContigs):
        if sc % 30 == 0:
            print(int(sc / len(listOfSuperContigs) * 1000) / 10, "% of exporting done")

        for c, contig in enumerate(supercontig):
            f.write("S\t" + names[contig] + "-" + str(copiesUsed[contig]) + "\t")
            if gfaFile != "":
                f.write(get_contig_GFA(gfaFile, names[contig]) + "\n")
            else:
                f.write("*\n")
            #print(supercontig)
            if c == 0:
                addresses[sc * 2] = names[contig] + "-" + str(copiesUsed[contig])
            if c == len(supercontig) - 1:
                addresses[sc * 2 + 1] = names[contig] + "-" + str(copiesUsed[contig])

            if c > 0:
                o1,o2 = -1,-1
                
                f.write("L\t"+ names[supercontig[c - 1]]+ "-"+ str(copiesUsed[supercontig[c - 1]] - 1))
                
                if contig in [int(i/2) for i in originalLinks[supercontig[c-1]*2]] :                    
                    f.write("\t-\t")
                    o1 = 0
                elif contig in [int(i/2) for i in originalLinks[supercontig[c-1]*2+1]]:
                    f.write("\t+\t")
                    o1 = 1
                else :
                    print('Problem while exporting, previously non-existing links seem to have been made up')
                    
                f.write(names[contig] + "-"+ str(copiesUsed[contig]))
                
                if supercontig[c-1] in [int(i/2) for i in originalLinks[contig*2]] :                    
                    f.write("\t+\t")
                    o2 = 0
                elif supercontig[c-1] in [int(i/2) for i in originalLinks[contig*2+1]]:
                    f.write("\t-\t")
                    o2 = 1
                else :
                    print('Problem while exporting, previously non-existing links seem to have been made up')
                
                if originalLinksCIGAR == None :
                    f.write("*\n")
                else :
                    f.write(originalLinksCIGAR[2*supercontig[c - 1]+o1][originalLinks[2*supercontig[c-1]+o1].index(2*contig+o2)]+'\n')

            copiesUsed[contig] += 1

    #then write in file the links between the ends of supercontigs
    #this part actually makes mistakes when two supercontigs are linked more than one time
    for endOfSuperContig in range(len(links)):
        for l in links[endOfSuperContig]:
            if endOfSuperContig < l:  # to write each link only once
                if l == endOfSuperContig + 1 and l%2 == 1: #if a supercontig loops on itself
                    f.write("L\t"+ addresses[endOfSuperContig]+'\t+\t'+ addresses[l] + '\t+\t')
                    if originalLinksCIGAR == None :
                        f.write('*\n')
                    else :       
                        if 2*listOfSuperContigs[int(l/2)][-(l%2)] in originalLinks[2*listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]]\
                            or 2*listOfSuperContigs[int(l/2)][-(l%2)] in originalLinks[2*listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]+1]:
                                o1 = 0
                        else :
                            o1 = 1
                            
                        if 2*listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)] in originalLinks[2*listOfSuperContigs[int(l/2)][-(l%2)]+o1] :
                            o2 = 0
                        else :
                            o2 = 1
                        endOfContig1 = 2*listOfSuperContigs[int(l/2)][-(l%2)]+o1
                        endOfContig2 = 2*listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]+o2
                        
                        f.write(originalLinksCIGAR[endOfContig1][originalLinks[endOfContig1].index(endOfContig2)] +'\n')
    
                else :
                    o1, o2 = -1, -1
                    
                    f.write("L\t"+ addresses[endOfSuperContig])
                    
                    if listOfSuperContigs[int(l/2)][-(l%2)] in [int(i/2) for i in originalLinks[listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]*2]] :
                        f.write("\t-\t")
                        o1 = 0
                    elif listOfSuperContigs[int(l/2)][-(l%2)] in [int(i/2) for i in originalLinks[listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]*2+1]] :
                        f.write("\t+\t")
                        o1 = 1
                    else :
                        print('Problem while exporting, previously non-existing links seem to have been made up')
                    
                        
                    f.write(addresses[l])
    
                    if listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)] in [int(i/2) for i in originalLinks[listOfSuperContigs[int(l/2)][-(l%2)]*2]] :
                        f.write("\t+\t")
                        o2 = 0
                    
                    elif listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)] in [int(i/2) for i in originalLinks[listOfSuperContigs[int(l/2)][-(l%2)]*2+1]] :
                        f.write("\t-\t")
                        o2 = 1
                    else :
                        print('Problem while exporting, previously non-existing links seem to have been made up')
                       
                    
                    if originalLinksCIGAR == None :
                        f.write('*\n')
                    else :
                        endOfContig1 = 2*listOfSuperContigs[int(l/2)][-(l%2)]+o2
                        endOfContig2 = 2*listOfSuperContigs[int(endOfSuperContig/2)][-(endOfSuperContig%2)]+o1
                        try :
                            f.write(originalLinksCIGAR[endOfContig1][originalLinks[endOfContig1].index(endOfContig2)] +'\n')
                        except :
                            f.write('*\n')
    




