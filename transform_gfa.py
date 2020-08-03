# -*- coding: latin-1 -*-
"""
This file is for handling the gfa : loading it into python lists, transforming it into fasta, reading it...
"""
import time
import sys

from segment import Segment


def gfa_to_fasta(gfaFilename="data/Assembly.gfa", fastaFilename="data/Assembly.fasta"):

    t1 = time.time()

    gfa_file = open(gfaFilename, "r")
    fasta_file = open(fastaFilename, "w")

    # line count
    i = 1
    seq_i = 0

    seqs = []

    for line in gfa_file.readlines():
        line = line.strip().split()
        if "S" in line[0]:
            if len(line) >= 3:
                fasta_file.write(">{0}\n{1}\n".format(line[1], line[2]))
                seqs.append(line[2])
                seq_i = seq_i + 1
            else:
                print("Wrong format in line {0}: expected three fields.".format(i))
                sys.exit(1)
        i = i + 1

    gfa_file.close()
    fasta_file.close()

    print("Processed {0} sequences.".format(seq_i))
    print(time.time() - t1)

    return seqs


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
            names[s.names[0]] = index #that is the (strange) way of adding a key to a dict in python
            index += 1
            

    print('Loading links')
    gfa_read = open(file, "r")
        
    for line in gfa_read:
        if line[0] == "L":

            l = line.strip('\n').split("\t")
            segments[names[l[1]]].add_link_from_GFA(line, names, segments, 0)
            segments[names[l[3]]].add_link_from_GFA(line, names, segments, 1)

    gfa_read.close()

    return segments, names

#print_short is useful to read a sequence file, shortening the sequences
def print_short():

    gfa_read = open("results/A_Vaga_PacBio/A_Vaga_finished2.gfa")
    r = gfa_read.read()

    bases = ["A", "C", "G", "T"]
    count = 0
    s = ""

    for i in range(1000000):
        if r[i] in bases:
            if count % 5000 == 0:
                s += r[i]
        else:
            s += r[i]
        count += 1
    print(s)
    return 0


# a function to test that links are, as they should, represented once at each of their extremities
def check_segments(listOfSegments):
    
    for segment in listOfSegments :
        for endOfSegment in range(2) :
            for n, neighbor in enumerate(segment.links[endOfSegment]) :
                if not segment in neighbor.links[segment.otherEndOfLinks[n]] :
                    print("Problem in links, a one-end link going from: ", segment.names, ' to ', neighbor.names)
                    return False
    return True

#function if you want to strip the suffix containing the copiesNumber
def strip_copiesNumber(gfaFileIn, gfaFileOut):
    with open(gfaFileIn, 'r') as f :
        with open(gfaFileOut, 'w') as fo:
            for line in f :
                
                l = line.split('\t')
                if 'S' in l[0] :
                    ll = l[1].split('-')
                    l[1] = ll[0]
                    
                if 'L' in l[0] :
                    ll = l[1].split('-')
                    l[1] = ll[0]
                    
                    ll = l[3].split('-')
                    l[3] = ll[0]
                    
                fo.write('\t'.join(l))
                
       
#strip_copiesNumber('Arabidopsis/Arabidopsis_hybrid/simplified_graph.gfa', 'Arabidopsis/Arabidopsis_hybrid/simplified_graph2.gfa')