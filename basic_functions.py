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
import re #to find all numbers in a mixed number/letters string (such as 31M1D4M)

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
    content = [x.strip("\n").split('\t') for x in content[1:]]
    # 1: contig_id, 2: length, 3: n_frags, 4:cumul_length
    content = [[x[0], int(x[1]), int(x[2]), int(x[3])] for x in content]

    return content


def interactionMatrix(hiccontactsfile, fragmentList, names, segments, header=True):  # the header refers to the hiccontactsfile

    print('Building the interaction matrix')
    t = time.time()
    # create interaction matrix of contig vs contig
    # 1 -> [1...N] N contigs
    # ...
    # N -> [1...N]
    interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))


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
        if contig1 in names and contig2 in names :
            index1 = names[contig1]
            index2 = names[contig2]
            
            if index1 != index2 :
                # add contacts to interaction matrix
                interactionMatrix[index1,index2] += contact[2]
                interactionMatrix[index2,index1] += contact[2]
                
                #adds the HiC coverage to the right contigs
                segments[index1].HiCcoverage += contact[2]
                segments[index2].HiCcoverage += contact[2]

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

#input : contig ID, gfa file and contigOffset, the position of the contig in the GFA file
#output : sequence
def get_contig_GFA(gfaFile, contig, contigOffset):
       
    with open(gfaFile) as f:

        f.seek(contigOffset)
        line = f.readline()         
        sline = line.strip('\n').split('\t')
        if len(sline) == 3 and sline[0] == 'S' and (contig in sline[1]) :
                return sline[2]
        elif len(sline) > 3  and sline[0] == 'S' and (contig in sline[1]) :
            return sline[2]+'\t'+sline[3]
        else :
            print('ERROR : Problem in the offset file, not pointing to the right lines')

    return "In get_contig : the contig you are seeking is not in the gfa file"

# Input :
#   offset file is for speeding up exportation
#   merge_adjacent_contig is to produce a GFA with contigs merged
def export_to_GFA(listOfSegments, gfaFile="", exportFile="results/newAssembly.gfa", offsetsFile = "", merge_adjacent_contigs = False): 
    
    #compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig
    noOffsets = False
    print('Offsets  : ', offsetsFile)
    if offsetsFile == "" :
        noOffsets = True
        offsetsFile = gfaFile.strip('.gfa') + '_offsets.pickle'
        
    if gfaFile != "" and noOffsets:
        #print("coucou")
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile :
            for line in gfafile:
                sline = line.strip('\n').split('\t')
                if sline[0] == 'S' :
                    #print('In export_to_GFA : exporting ', sline[1])
                    line_offset[sline[1]] = offset #adds pair sline[1]:offset to the dict
                    
                offset += len(line)
            
        with open(offsetsFile, 'wb') as o:
            pickle.dump(line_offset, o)
    
    if gfaFile != '' :
        with open(offsetsFile, 'rb') as o:
            line_offset = pickle.load(o)
        
        #print(line_offset)
 
    print('Line_offsets loaded, launching proper writing of the new GFA')
    #Now that the preliminary work is done, start writing the new gfa file    

    f = open(exportFile, "w")
    
    #compute the copiesnumber
    compute_copiesNumber(listOfSegments)

    #write the sequences and the links within the supercontigs
    t = time.time()
    
    if merge_adjacent_contigs == False :
        for s, segment in enumerate(listOfSegments):
            if  time.time() > t+1 :
                t = time.time()
                print(int(s / len(listOfSegments) * 1000) / 10, "% of sequences written", end = '\r')
    
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
    
        print('Done exporting sequences, just a little more time...')
        #then write in the gfa file the links between the ends of supercontigs
    
        for s, segment in enumerate(listOfSegments):
            
            if  time.time() > t+1 :
                t = time.time()
                print(int(s / len(listOfSegments) * 1000) / 10, "% of links written", end = '\r')
                
            for endOfSegment in range(2):
                for l, neighbor in enumerate(segment.links[endOfSegment]):
                    
                    if segment.ID <= neighbor.ID : #that is to ensure each link is written only once
                    
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

    # in the case the user prefers having merged contigs as an output
    else : #if merge_adjacent_contigs == True
        
        for s, segment in enumerate(listOfSegments):
            if  time.time() > t+1 :
                t = time.time()
                print(int(s / len(listOfSegments) * 1000) / 10, "% of sequences written", end = '\r')
            
            f.write("S\t" + segment.full_name() + "\t")
            if gfaFile != "":
                
                sequence = ''
                for c, contig in enumerate(segment.names) :
                    s = get_contig_GFA(gfaFile, contig, line_offset[contig])
                    if segment.orientations[c] == 0 :
                        s = s[::-1]
                    if c > 0 :
                        CIGARlength = np.sum([int(i) for i in re.findall(r'\d+', segment.insideCIGARs[c-1])])
                        s = s[CIGARlength:]
                    sequence += s
                f.write(sequence + "\n")
                
            else:
                f.write("*\n")
                
            for endOfSegment in range(2) :
                for n, neighbor in enumerate(segment.links[endOfSegment]):
                    if segment.ID < neighbor.ID : #to write each link just one
                        orientation1, orientation2 = '+', '+'
                        if endOfSegment == 0 :
                            orientation1 = '-'
                        if segment.otherEndOfLinks[endOfSegment][n] == 1 :
                            orientation2 = '-'
                            
                        f.write("L\t"+segment.full_name()+'\t'+orientation1+'\t'+neighbor.full_name()+\
                                '\t'+orientation2+'\t'+ segment.CIGARs[endOfSegment][n]+'\n')




