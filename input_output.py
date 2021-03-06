#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""
import numpy as np
from scipy import sparse #to handle interactionMatrix, which should be sparse
import time #to inform the user on what the programm is doing on a regular basis
import os.path #to check the existence of files
import pickle #for writing files and reading them
import re #to find all numbers in a mixed number/letters string (such as 31M1D4M), to split on several characters (<> in longReads_interactionMatrix)

from segment import Segment
from segment import compute_copiesNumber
from segment import delete_links_present_twice

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
    unknowncontacts = 0
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
                
        else :
            unknowncontacts += 1

        n += 1
        
    if unknowncontacts != 0 :
        print('There are ', unknowncontacts, ' out of ', n, ' contacts I did not manage to map : you may want to check if the names of the contigs are consistent throughout your files')
        
    return interactionMatrix

#input : a gaf file (outputted by graphaligner)
#output : an interaction matrix with two values : 100 if the contigs are next to each other in some reads, 0 elsewhise
def longReads_interactionsMatrix(gafFile, names, segments, similarity_threshold = 0, whole_mapping = False):
     
    print('Building interaction matrix from the gaf file')
    f = open(gafFile, 'r')
    
    interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))
    repeats = [0]*len(segments)
    
    allLinks = set()
    noNames = set() #list of names that are not found in names
   
    #these two values are in case some contigs of the gafs are just subnames of names
    namesPlus = names.copy()
    contigsPlus = {}
    for i in names :
        contigsPlus[i] = i
    
    for line in f :
        
        ls = line.split('\t')
        if ls[5].count('>') + ls[5].count('<') > 1 :
            
        
            if (not 'id:f' in ls[-2]) or (float(ls[-2].split(':')[-1]) > similarity_threshold) :
                
                if not whole_mapping or (float(ls[2]) == 0 and float(ls[1]) == float(ls[3])) :
                    
                    contigs = re.split('[><]' , ls[5])
                    orientations = "".join(re.findall("[<>]", ls[5]))
                    del contigs[0] #because the first element is always ''

                    for c1 in range(len(contigs)-1) :
                        for c2 in range(c1+1, len(contigs)):
                            
                            if contigs[c1] not in noNames and contigs[c2] not in noNames :
                            
                                #try to look for the name if it is not in names
                                if contigs[c1] not in namesPlus :
                                    found = False
                                    for i in names.keys() :
                                        if contigs[c1] in i :
                                            namesPlus[contigs[c1]] = names[i]
                                            contigsPlus[contigs[c1]] = i
                                            found = True
                                    if not found :
                                        noNames.add(contigs[c1]) #if you do not find it once, the cause is lost
                                            
                                if contigs[c2] not in namesPlus :
                                    found = False
                                    for i in names.keys() :
                                        if contigs[c2] in i :
                                            namesPlus[contigs[c2]] = names[i]
                                            contigsPlus[contigs[c2]] = i
                                            found = True
                                    if not found :
                                        noNames.add(contigs[c2])
                                
                                #if the two contigs exist, go ahead
                                if contigs[c1] in namesPlus and contigs[c2] in namesPlus :

                                    if namesPlus[contigs[c1]] != namesPlus[contigs[c2]] :
                                        interactionMatrix[namesPlus[contigs[c1]], namesPlus[contigs[c2]]] += 1
                                        if c2 == c1 +1 :
                                            if contigsPlus[contigs[c1]] == contigs[c1] and contigsPlus[contigs[c2]] == contigs[c2] :

                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                                
                                            elif contigsPlus[contigs[c1]] == contigs[c1] :
                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], False))
                                            elif  contigsPlus[contigs[c2]] == contigs[c2] :
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                            else :
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], False))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], False))
                                    else :
                                        repeats[namesPlus[contigs[c1]]] = max(contigs.count(contigs[c1]), repeats[namesPlus[contigs[c1]]])
                    
    return interactionMatrix, repeats, allLinks

def load_interactionMatrix(file, listOfSegments, names) :
    f = open(file, 'rb')
    interactionMatrix  = pickle.load(f)
    for segment in listOfSegments :
        for contig in segment.names :
            segment.HiCcoverage += np.sum(interactionMatrix[names[contig]])
    
    return interactionMatrix

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
#output : sequence, and if it is present, the sequencing depth of the contig
def get_contig_GFA(gfaFile, contig, contigOffset):
       
    with open(gfaFile) as f:

        f.seek(contigOffset)
        line = f.readline()         
        sline = line.strip('\n').split('\t')
        if len(sline) == 3 and sline[0] == 'S' and (contig in sline[1]) :
                return sline[2], ''
            
        elif len(sline) > 3 and sline[0] == 'S' and (contig in sline[1]):
            for f in sline[3:] :
                if 'dp' in f or 'DP' in f :
                    return sline[2], f
            return sline[2], ''

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
            
        # with open(offsetsFile, 'wb') as o:
        #     pickle.dump(line_offset, o)
    
    # if gfaFile != '' :
    #     with open(offsetsFile, 'rb') as o:
    #         line_offset = pickle.load(o)
        
        #print(line_offset)
 
    print('Line_offsets computed, launching proper writing of the new GFA')
    #Now that the preliminary work is done, start writing the new gfa file    

    # for s in listOfSegments :


    f = open(exportFile, "w")
    
    #compute the copiesnumber
    copies = compute_copiesNumber(listOfSegments)

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
                    sequence, depth = get_contig_GFA(gfaFile, contig, line_offset[contig])
                    if depth == '':
                        f.write(sequence + "\n")
                    else :
                        newdepth = str(float(depth.split(':')[-1])/copies[contig])
                        f.write(sequence + '\tDP:f:'+newdepth + '\n')
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
            
            fullDepth = 0
            
            if gfaFile != "":
                
                sequence = ''
                for c, contig in enumerate(segment.names) :
                    s, depth = get_contig_GFA(gfaFile, contig, line_offset[contig])
                    if segment.orientations[c] == 0 :
                        s = s[::-1]
                    if c > 0 :
                        CIGARlength = np.sum([int(i) for i in re.findall(r'\d+', segment.insideCIGARs[c-1])])
                        
                        s = s[CIGARlength:]
                    if depth != '' :
                        fullDepth += ( float(depth.split(':')[-1])/copies[contig] ) * len(s)
                    
                    sequence += s
                    
                if fullDepth == 0:
                    f.write(sequence + "\n")
                else :
                    newdepth = str(fullDepth/len(sequence))
                    f.write(sequence + '\tDP:f:'+newdepth + '\n')
                
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
                                
def export_to_fasta(listOfSegments, gfaFile, exportFile="results/newAssembly.fasta"): 
    
    #compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig

    t = 0
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
            
 
    print('Line_offsets computed, launching writing of the fasta')
    #Now that the preliminary work is done, start writing the new fasta file    

    f = open(exportFile, "w")
    
    #compute the copiesnumber
    copies = compute_copiesNumber(listOfSegments)


    # now sort the segments by length, to output at the beginning of the files the longests fragments
    listOfSegments.sort(key = lambda x : x.length, reverse = True)
    
    
    #Finally, write the sequences
    for s, segment in enumerate(listOfSegments):
        
        if  time.time() > t+1 :
            t = time.time()
            print(int(s / len(listOfSegments) * 1000) / 10, "% of sequences written", end = '\r')
        
        f.write(">supercontig_" + str(s+1) + "\n")
        
        fullDepth = 0
        

        
        sequence = ''
        for c, contig in enumerate(segment.names) :
            s, depth = get_contig_GFA(gfaFile, contig, line_offset[contig])
            if segment.orientations[c] == 0 :
                s = s[::-1]
            if c > 0 :
                CIGARlength = np.sum([int(i) for i in re.findall(r'\d+', segment.insideCIGARs[c-1])])
                
                s = s[CIGARlength:]
            if depth != '' :
                fullDepth += ( float(depth.split(':')[-1])/copies[contig] ) * len(s)
            
            sequence += s
            
        f.write(sequence + "\n")


# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
# Also returns the list of the contig's names
def load_gfa(file):

    print('Loading contigs')
    gfa_read = open(file, "r")

    segments = []
    
    index = 0
    names = {}
    for line in gfa_read:
        if line[0] == "S":
            l = line.strip('\n').split("\t")
            s = Segment([l[1]], [1], [len(l[2])])
            segments.append(s)
            names[s.names[0]] = index #that is the (strange) way of adding a new key to a dict in python
            index += 1
            

    print('Loading links')
    gfa_read = open(file, "r")
        
    for line in gfa_read:
        if line[0] == "L":

            l = line.strip('\n').split("\t")
            segments[names[l[1]]].add_link_from_GFA(line, names, segments, 0)
            segments[names[l[3]]].add_link_from_GFA(line, names, segments, 1)

    gfa_read.close()
    
    delete_links_present_twice(segments)

    return segments, names

