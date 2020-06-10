#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""
import pandas as pd
import numpy as np
from scipy import sparse #to handle interactionMatrix, which should be sparse
import time #to inform the user on what the programm is doing on a regular basis
import os.path #to check the existence of files
import pickle #for writing files and reading them

from segment import Segment
from segment import compute_copiesNumber

# Read fragments list file
# Input :
#   file : fragments_list.txt
# Output :
#   content : list with contig_id, fragment start, fragment end, fragment length
def read_fragment_list(file, header=True):

    with open(file) as f:
        content = f.readlines()

    if header:
        del content[0]

    # parsing
    # 1: contig_id, 2: fragment_start, 3: fragment_end, 4: fragment_length
    content = [x.strip("\n").split("\t") for x in content]
    content = [[x[1], int(x[2]), int(x[3]), int(x[4])] for x in content]

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

    print('Building the interaction matrix')
    t = time.time()
    # create interaction matrix of contig vs contig
    # 1 -> [1...N] N contigs
    # ...
    # N -> [1...N]
    interactionMatrix = sparse.dok_matrix((len(names), len(names)))


    with open(hiccontactsfile) as f:
        inFile = f.readlines()

    if header:
        del inFile[0]
    
    n = 0
    for line in inFile:
        
        if time.time()-t > 2 :
            t = time.time()
            print('Built '+str(int(n/len(inFile)*100))+'%', end='\r')

        line = line.strip("\n").split("\t")

        # frag1, frag2, contacts
        contact = [int(line[0]), int(line[1]), int(line[2])]

        # search for contig name corresponding to fragment id
        contig1 = fragmentList[contact[0]][0]
        contig2 = fragmentList[contact[1]][0]

        # search for the index of the contigs in names 
        index1 = names[contig1]
        index2 = names[contig2]
        
        # add contacts to interaction matrix
        interactionMatrix[index1,index2] += contact[2]
        interactionMatrix[index1,index2] += contact[2]

        n += 1
    return interactionMatrix


def export_to_csv(l, file):
    df = pd.DataFrame(l)
    df.to_csv(file)


def import_from_csv(file):

    df = pd.read_csv(file)
    l = df.values.tolist()

    newl = [[x for x in i if not pd.isnull(x)] for i in l]
    return [x[1:] for x in newl]  # we discard the header line

#input : contig ID and fasta file
#output : sequence
def get_contig_FASTA(fastaFile, contig, firstline=0):

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

#input : contig ID and gfa file
#output : sequence
def get_contig_GFA(gfaFile, contig, contigOffset):
    
    with open(gfaFile) as f:

        f.seek(contigOffset)
        line = f.readline()         
        sline = line.strip('\n').split('\t')
        if len(sline) >= 3 and sline[0] == 'S' and (contig in sline[1]) :
                return sline[2]
        else :
            print('ERROR : Problem in the offset file, not pointing to the right lines')

    return "In get_contig : the contig you are seeking is not in the gfa file"

def export_to_GFA(listOfSegments, gfaFile="", exportFile="results/newAssembly.gfa", offsetsFile = ""): #offset file is for speeding up exportation : 

    #compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig
    if offsetsFile == "" :
        offsetsFile = gfaFile.strip('.gfa') + '_offsets.pickle'
        
    if gfaFile != "" and not os.path.exists(offsetsFile) :
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile :
            for line in gfafile:
                sline = line.strip('\n').split('\t')
                if sline[0] == 'S' :
                    line_offset[sline[1]] = offset #adds pair sline[1]:offset to the dict
                    
                offset += len(line)
            
        with open(offsetsFile, 'wb') as o:
            pickle.dump(line_offset, o)
    
    if gfaFile != '' :
        with open(offsetsFile, 'rb') as o:
            line_offset = pickle.load(o)
        
        print(line_offset)    
        

    f = open(exportFile, "w")
    
    #compute the copiesnumber
    cn = compute_copiesNumber(listOfSegments)

    #write the sequences and the links within the supercontigs
    for s, segment in enumerate(listOfSegments):
        if s % 30 == 0:
            print(int(s / len(listOfSegments) * 1000) / 10, "% of exporting done")

        for c, contig in enumerate(segment.names):
            
            f.write("S\t" + contig + "-" + str(segment.copiesnumber[c]) + "\t")
            if gfaFile != "":
                f.write(get_contig_GFA(gfaFile, contig, line_offset[contig]) + "\n")
            else:
                f.write("*\n")

            if c > 0:
                
                f.write("L\t"+ segment.names[c-1]+ "-"+ str(segment.copiesnumber[c-1]))
                
                if segment.orientations[c-1] == 1 :                    
                    f.write("\t+\t")

                elif segment.orientations[c-1] == 0:
                    f.write("\t-\t")
             
                f.write(contig + "-"+ str(segment.copiesnumber[c]))
                
                if segment.orientations[c] == 1 :                    
                    f.write("\t+\t")

                elif segment.orientations[c] == 0:
                    f.write("\t-\t")
                    
                f.write(segment.insideCIGARs[c-1]+'\n')

    #then write in the gfa file the links between the ends of supercontigs

    for s, segment in enumerate(listOfSegments):
        for endOfSegment in range(2):
            for l, neighbor in enumerate(segment.links[endOfSegment]):
                
                if segment.hash() <= neighbor.hash() : #that is to ensure each link is written only once
                
                    endOfNeighbor = segment.otherEndOfLinks[endOfSegment][l]
                    orientation1, orientation2 = '-', '-'
                    
                    if segment.orientations[-endOfSegment] == endOfSegment :
                        orientation1 = '+'
                        
                    if neighbor.orientations[-endOfNeighbor] != endOfNeighbor :
                        orientation2 = '+'
                        
                    f.write("L\t"+segment.names[-endOfSegment] +"-"+ str(segment.copiesnumber[-endOfSegment]) + '\t' \
                            + orientation1 + '\t' +\
                                neighbor.names[-endOfNeighbor] +"-"+ str(neighbor.copiesnumber[-endOfNeighbor])+'\t'\
                            +orientation2+'\t'+segment.CIGARs[endOfSegment][l]+'\n')
    




