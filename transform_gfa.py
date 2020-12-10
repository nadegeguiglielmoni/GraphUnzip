# -*- coding: latin-1 -*-
"""
This file is for handling the gfa : loading it into python lists, transforming it into fasta, reading it...
"""
import time
import sys

from segment import Segment


def gfa_to_fasta(gfaFilename, fastaFilename = ''):

    if fastaFilename == '' :
        fastaFilename = gfaFilename.strip('gfa')+'fasta'
        
    t1 = time.time()

    gfa_file = open(gfaFilename, "r")
    fasta_file = open(fastaFilename, "w")

    # line count
    i = 1
    seq_i = 0

    seqs = []

    for line in gfa_file.readlines():
        line = line.strip().split()
        if len(line) > 3 :
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
    print(time.time() - t1, 's')

    return seqs


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

#function if you want to strip the suffix containing the copiesNumber. Careful though, if contigs are duplicated they will be merged back since they will have the same name
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
#gfa_to_fasta('Arabidopsis/Arabidopsis_hybrid/assembly_graph.gfa')
#gfa_to_fasta("Escherichia_Coli/1a1k/assemblyGraph_k63_noOverlaps.gfa")
#gfa_to_fasta('bacteria_mix/SPAdes_output/assembly_graph_after_simplification.gfa')
#gfa_to_fasta('data_A_vaga_HiFi/Flye/A_vaga_12_chr.gfa')