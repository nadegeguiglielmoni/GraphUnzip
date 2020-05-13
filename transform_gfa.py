# -*- coding: latin-1 -*-
"""
This file is for handling the gfa : loading it into python lists, transforming it into fasta, reading it...
"""
import time
import sys


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

    names = []
    for line in gfa_read:
        if line[0] == "S":
            l = line.split("\t")
            names.append(l[1])

    links = [
        [] for i in range(len(names) * 2)
    ]  # a list of links at one end of a contig

    for line in gfa_read:
        if line[0] == "L":

            if len(line) < 5:
                print(
                    "Wrong format: expected at least 5 fields in line:\n{0}\n".format(
                        line
                    )
                )

            l = line.split("\t")
            contig1 = l[1]
            contig2 = l[3]

            contig1index = names.index(contig1)
            contig2index = names.index(contig2)

            if l[2] == "+" and l[4] == "+":
                links[contig1index * 2 + 1] += [contig2index * 2]
                links[contig2index * 2] += [contig1index * 2 + 1]

            elif l[2] == "+" and l[4] == "-":
                links[contig1index * 2 + 1] += [contig2index * 2 + 1]
                links[contig2index * 2 + 1] += [contig1index * 2 + 1]

            elif l[2] == "-" and l[4] == "-":
                links[contig1index * 2] += [contig2index * 2 + 1]
                links[contig2index * 2 + 1] += [contig1index * 2]

            elif l[2] == "-" and l[4] == "+":
                links[contig1index * 2] += [contig2index * 2]
                links[contig2index * 2] += [contig1index * 2]
            else:
                print("There seems to be a problem in the gfa file.")

    gfa_read.close()

    return links, names


def print_short():

    gfa_read = open("data/Assembly.gfa")
    r = gfa_read.read()

    bases = ["A", "C", "G", "T"]
    count = 0
    s = ""

    for i in range(len(r)):
        if r[i] in bases:
            if count % 5000000 == 0:
                s += r[i]
        else:
            s += r[i]
        count += 1
    print(s)
    return 0


# print_short()


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
