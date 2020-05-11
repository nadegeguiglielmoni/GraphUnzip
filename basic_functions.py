#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

@author: zaltabar

File basically dedicated to small functions involving reading and writing files
"""
import pandas as pd
import numpy as np

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


def export_to_csv(l, file):
    df = pd.DataFrame(l)
    df.to_csv(file)


def import_from_csv(file):

    df = pd.read_csv(file)
    l = df.values.tolist()

    newl = [[x for x in i if not np.isnan(x)] for i in l]
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


def export_to_GFA(
    links,
    listOfSuperContigs,
    copiesnumber,
    names=None,
    fastaFile="",
    exportFile="results/newAssembly.gfa",
):

    if names == None:
        names = [2 * i for i in range(len(links) / 2)]
    #  copyfile("results/sequencesWithoutLinks.gfa", "results/newAssembly.gfa")
    f = open(exportFile, "w")
    addresses = ["" for i in range(len(listOfSuperContigs) * 2)]
    copiesUsed = [0 for i in copiesnumber]

    for sc, supercontig in enumerate(listOfSuperContigs):
        if sc % 30 == 0:
            print(int(sc / len(listOfSuperContigs) * 1000) / 10, "% of exporting done")
        for c, contig in enumerate(supercontig):
            f.write("S\t" + names[contig] + "-" + str(copiesUsed[contig]) + "\t")
            if fastaFile != "":
                f.write(get_contig(fastaFile, contig, contig * 2 - 1) + "\n")
            else:
                f.write("*\n")
            # print(supercontig)
            if c == 0:
                addresses[sc * 2] = names[contig] + "-" + str(copiesUsed[contig])
            if c == len(supercontig) - 1:
                addresses[sc * 2 + 1] = names[contig] + "-" + str(copiesUsed[contig])

            if c > 0:
                # /!\ the + orientation of both contigs next line is arbitrary, do NOT trust it to build an actual genome
                f.write(
                    "L\t"
                    + names[supercontig[c - 1]]
                    + "-"
                    + str(copiesUsed[supercontig[c - 1]] - 1)
                    + "\t+\t"
                    + names[contig]
                    + "-"
                    + str(copiesUsed[contig])
                    + "\t+\t*\n"
                )

            copiesUsed[contig] += 1

    for endOfSuperContig in range(len(links)):
        for l in links[endOfSuperContig]:
            if endOfSuperContig < l:  # to write each link only once
                if endOfSuperContig % 2 == 1 and l % 2 == 0:
                    f.write(
                        "L\t"
                        + addresses[endOfSuperContig]
                        + "\t+\t"
                        + addresses[l]
                        + "\t+\t*\n"
                    )
                elif endOfSuperContig % 2 == 0 and l % 2 == 0:
                    f.write(
                        "L\t"
                        + addresses[endOfSuperContig]
                        + "\t-\t"
                        + addresses[l]
                        + "\t+\t*\n"
                    )
                elif endOfSuperContig % 2 == 0 and l % 2 == 1:
                    f.write(
                        "L\t"
                        + addresses[endOfSuperContig]
                        + "\t-\t"
                        + addresses[l]
                        + "\t-\t*\n"
                    )
                elif endOfSuperContig % 2 == 1 and l % 2 == 1:
                    f.write(
                        "L\t"
                        + addresses[endOfSuperContig]
                        + "\t+\t"
                        + addresses[l]
                        + "\t-\t*\n"
                    )

