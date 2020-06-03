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

    gfa_read = open(file, "r")

    segments = []
    
    for line in gfa_read:
        if line[0] == "S":
            l = line.strip('\n').split("\t")
            s = Segment(len(segments),[len(segments)], [l[1]], [1], [len(l[2])])
            segments.append(s)

    
    gfa_read = open(file, "r")
        
    names = [i.names[0] for i in segments]
    print(names)
    for line in gfa_read:
        if line[0] == "L":

            l = line.strip('\n').split("\t")
            
            segments[names.index(l[1])].add_link_from_GFA(line, names, segments)
            segments[names.index(l[3])].add_link_from_GFA(line, names, segments)

    gfa_read.close()

    return segments


def print_short():

    gfa_read = open("results/A_Vaga_PacBio/A_Vaga_finished2.gfa")
    r = gfa_read.read()

    bases = ["A", "C", "G", "T"]
    count = 0
    s = ""

    for i in range(1000000):
        if r[i] in bases:
            if count % 5000000 == 0:
                s += r[i]
        else:
            s += r[i]
        count += 1
    print(s)
    return 0

def print_short_gfa():
    
    with open("results/A_Vaga_PacBio/A_Vaga_finished2.gfa", 'r') as f :
        for line in f :
            print(line.split('\t')[0], line.split('\t')[1])

# a function to test that load_gfa worked ok
def check_links(links):
    for i, noeud in enumerate(links):
        for j, neighbor in enumerate(noeud):
            found = False
            for k in links[neighbor]:
                if k == i:
                    found = True
            if not found and neighbor != -1:
                print("Problem in links, a one-end link : ", i, neighbor)
                return False
    return True

#gfa_to_fasta('data_A_Vaga_Illumina/assemblyGraph_k201.gfa', 'data_A_Vaga_Illumina/assemblyGraph_k201.fasta')