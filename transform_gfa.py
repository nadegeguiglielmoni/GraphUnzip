# -*- coding: latin-1 -*-
"""
This file is for handling the gfa : loading it into python lists, transforming it into fasta, reading it...
"""
import time

# from analyse_HiC import lookForStrings
# print(r[:100000])


def gfa_to_fasta():

    gfa_read = open("data/Assembly.gfa")
    r = gfa_read.read()

    seqs = []
    current = ""
    count = 0
    bases = ["A", "C", "G", "T"]

    t1 = time.time()

    for i in range(len(r)):
        count += 1
        if current == "":
            if (
                r[i] in bases
                and r[i + 1] in bases
                and r[i + 2] in bases
                and r[i + 3] in bases
            ):
                current = r[i]

        else:
            # if i%4 == 0 :
            if r[i] in bases:
                current += r[i]
            else:
                seqs += [current]
                current = ""
            # else :
            #   current += r[i]

        if count % 1000000 == 0:
            print(count / 1000000)

        if count == 123533473:
            # print(current)
            break

    fasta_file = open("data/Assembly.fasta", "w")

    for i in range(len(seqs)):
        name = ">sequence" + str(i)
        fasta_file.write(name + "\n")
        fasta_file.write(seqs[i] + "\n")

    print(len(seqs))
    print(len(r))
    print(time.time() - t1)

    return seqs


# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
def load_gfa(lenseqs):  # algo marche si les noms des séquences sont leur indice *2

    gfa_read = open("Assembly.gfa")
    r = gfa_read.read()

    count = 0
    digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    orientation = ["-", "+"]

    links = [[] for i in range(lenseqs * 2)]  # a list of links at one end of a contig
    number = ""
    n1 = -1  # integers that refer to the two contigs of interest
    n2 = -1
    n1orientation = -1  # refer to the orientation of the contig (0 if -, 1 if +)
    n2orientation = -1
    link = False

    for i in range(len(r)):
        count += 1
        if r[i] == "L":
            link = True

        if link == True:

            if r[i] in digits:
                number += r[i]

            if r[i] in orientation:

                n = int(number)
                number = ""

                if n1 == -1:
                    n1 = n
                    if r[i] == "-":
                        n1orientation = 0
                    else:
                        n1orientation = 1

                else:
                    n2 = n
                    if r[i] == "-":
                        n2orientation = 0
                    else:
                        n2orientation = 1

                    links[int(n1 + n1orientation)] += [int(n2 + 1 - n2orientation)]
                    links[int(n2 + 1 - n2orientation)] += [int(n1 + n1orientation)]

                    # resetting all the parameters
                    n1 = -1
                    n2 = -1
                    n1orientation = -1
                    n2orientation = -1
                    link = False

    return links


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


# a function to test that our gfa_to_python worked ok
def check_links(links):
    for i, noeud in enumerate(links):
        for j, neighbor in enumerate(noeud):
            found = False
            for k in links[neighbor]:
                if k == i:
                    found = True
            if not found:
                print("Problem : ", i, neighbor)


print_short()
