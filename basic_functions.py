#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""
import pandas as pd
import numpy as np

from segment import Segment

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
        [x[1], int(x[2]), int(x[3]), int(x[4])] for x in content
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


def interactionMatrix(hiccontactsfile, fragmentList, names, header=True):  # the header refers to the hiccontactsfile

    # create a full interaction matrix of contig vs contig
    # 1 -> [1...N] N contigs
    # ...
    # N -> [1...N]
    interactionMatrix = [[0 for i in range(len(names))] for j in range(len(names))]

    with open(hiccontactsfile) as f:
        inFile = f.readlines()

    if header:
        del inFile[0]

    for line in inFile:

        line = line.strip("\n").split("\t")

        # frag1, frag2, contacts
        contact = [int(line[0]), int(line[1]), int(line[2])]

        # search for contig name corresponding to fragment id
        contig1 = fragmentList[contact[0]][0]
        contig2 = fragmentList[contact[1]][0]

        # search for the index of the contigs in names 
        index1 = names.index(contig1)
        index2 = names.index(contig2)
        
        # add contacts to interaction matrix
        interactionMatrix[index1][index2] += contact[2]
        interactionMatrix[index1][index2] += contact[2]

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

def export_to_GFA(listOfSegments, copiesnumber, gfaFile="", exportFile="results/newAssembly.gfa"):

    f = open(exportFile, "w")
    copiesUsed = [0 for i in copiesnumber]

    #fist write the sequences and the links within the supercontigs
    for s, segment in enumerate(listOfSegments):
        if s % 30 == 0:
            print(int(s / len(listOfSegments) * 1000) / 10, "% of exporting done")

        for c, contig in enumerate(segment.names):
            
            f.write("S\t" + contig + "-" + str(copiesUsed[segment.listOfContigs[c]]) + "\t")
            if gfaFile != "":
                f.write(get_contig_GFA(gfaFile, contig) + "\n")
            else:
                f.write("*\n")

            if c > 0:
                
                f.write("L\t"+ segment.names[c-1]+ "-"+ str(copiesUsed[segment.listOfContigs[c-1]]))
                
                if segment.orientations[c-1] == 1 :                    
                    f.write("\t+\t")

                elif segment.orientations[c-1] == 0:
                    f.write("\t-\t")
             
                f.write(contig + "-"+ str(copiesUsed[segment.listOfContigs[c]]))
                
                if segment.orientations[c] == 1 :                    
                    f.write("\t+\t")

                elif segment.orientations[c] == 0:
                    f.write("\t-\t")
                    
                f.write(segment.insideCIGARs[c-1]+'\n')

            copiesUsed[segment.listOfContigs[c]] += 1

    #then write in the gfa file the links between the ends of supercontigs

    for s, segment in enumerate(listOfSegments):
        for endOfSegment in range(2):
            for l, neighbor in enumerate(segment.links[endOfSegment]):
                
                if segment.hash <= neighbor.hash : #that is to ensure each link is written only once
                
                    endOfNeighbor = segment.otherEndOfLinks[l]
                    orientation1, orientation2 = '-', '-'
                    
                    if segment.orientations[-endOfSegment] == endOfSegment :
                        orientation1 = '+'
                        
                    if neighbor.orientations[-endOfNeighbor] != endOfNeighbor :
                        orientation2 = '+'
                        
                    f.write("L\t"+segment.names[-endOfSegment] + '\t' + orientation1 + '\t' + neighbor.names[-endOfNeighbor]+'\t'\
                            +orientation2+'\t'+segment.CIGARs[endOfSegment][l]+'\n')
    




