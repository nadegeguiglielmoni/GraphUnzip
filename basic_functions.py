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
# input :
#   file : sparse matrix with 3 fields: frag1, frag2, contacts
# output :
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
def read_fragment_list(file):

    with open(file) as f:
        content = f.readlines()

    # parsing
    content = [x.strip("\n") for x in content[1:]]
    content = [x.split("\t") for x in content]
    content = [
        [int(x[1].strip("sequence")), int(x[2]), int(x[3]), int(x[4])] for x in content
    ]

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
                target = ">sequence" + str(contig)
                if target in line:
                    lookAtNextLine = True

            linenumber += 1
    return "In get_contig : the contig you are seeking is not in the fasta file"

