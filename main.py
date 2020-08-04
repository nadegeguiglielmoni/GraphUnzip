#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""

import basic_functions as bf
import analyse_HiC
from transform_gfa import load_gfa
from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities
from segment import check_if_all_links_are_sorted

from scipy import sparse
import numpy as np
import argparse
import os.path
import sys
import pickle #reading and writing files
import time

def parse_args():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gfa", required=True, help="""GFA file to phase""")
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output GFA [default: output.gfa]""",
    )
    parser.add_argument(
        "-A",
        "--accepted",
        required=False,
        default=0.20,
        help="""Threshold to accept links. [default: 0.20]""",
    )
    parser.add_argument(
        "-R",
        "--rejected",
        required=False,
        default=0.10,
        help="""Threshold to reject links. [default: 0.10]""",
    )
    parser.add_argument(
        "-s",
        "--steps",
        required=False,
        default=10,
        help="""Number of cycles get rid of bad links - duplicate contigs. [default: 10]""",
    )

    parser.add_argument(
        "-m", "--matrix", required=False, default="Empty", help="""Sparse contact map"""
    )
    parser.add_argument(
        "-F", "--fragments", required=False, default="Empty", help="""Fragments list"""
    )
    parser.add_argument(
        "-i",
        "--interactions",
        required=False,
        default="interactionMatrix.csv",
        help="""File with interactions [default: interactionMatrix.csv]""",
    )
    parser.add_argument(
        "--merge", required=False, default = "Empty",  help="""If you want the output to have all possible contigs merged (y/n) [default: n]"""
    )
    return parser.parse_args()

def main():

    args = parse_args()
    gfaFile = args.gfa
    outFile = args.output
    matrixFile = args.matrix
    fragmentsFile = args.fragments
    interactionFile = args.interactions
    stringenceReject = float(args.rejected)
    stringenceAccept = float(args.accepted)
    steps = int(args.steps)
    merge = args.merge
    
    t = time.time()

    if not os.path.exists(gfaFile):
        print("Error: could not find GFA file {0}.".format(gfaFile))
        sys.exit(1)

    # Loading the data
    print('Loading the GFA file')
    segments, names = load_gfa(gfaFile) #outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...

    if fragmentsFile is not "Empty" and matrixFile is not "Empty":
        if os.path.exists(fragmentsFile):
            fragmentList = bf.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix
            interactionMatrix = bf.interactionMatrix(matrixFile, fragmentList, names)

            # exporting it as to never have to do it again

            print('Exporting interaction matrix')
            with open(interactionFile, 'wb') as o:
                pickle.dump(interactionMatrix, o)

        else:
            print("Error: could not find fragments file {0}.".format(fragmentsFile))
            sys.exit(1)
    elif interactionFile is not "Empty" :
        print('Loading the interaction matrix')
        with open(interactionFile, 'rb') as o:
            interactionMatrix = pickle.load(o)
    else:
        if not os.path.exists(interactionFile):
            print("Error: you should provide either a processed interaction file, or the fragments list and the sparse contact map.")
            sys.exit(1)

    print('Everything loaded, moving on to solve_ambiguities')

    segments = solve_ambiguities(segments, interactionMatrix, stringenceReject, stringenceAccept, steps)

    # now exporting the output
    print('Now exporting')
    merge_adj = False
    if merge != "Empty" and merge != "n":
        merge_adj = True
    bf.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs = merge_adj)
    print('Finished in ', time.time()-t , " seconds")

# if __name__ == "__main__":
#     main()

gfaFile = "Arabidopsis/Arabidopsis_hybrid/simplified_graph.gfa"
#gfaFile = "Arabidopsis/Arabidopsis_hybrid/small2.gfa"
fragmentsFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/fragments_list.txt"
matrixFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/abs_fragments_contacts_weighted.txt"
interactionFile = "Arabidopsis/Arabidopsis_hybrid/unzip_out/interaction_matrix.pickle"
outFile = "Arabidopsis/Arabidopsis_hybrid/unzip_out/unzipped.gfa"

# gfaFile = "data_A_Vaga_PacBio/Assembly.gfa"
# #gfaFile = "Arabidopsis/Arabidopsis_hybrid/small2.gfa"
# fragmentsFile = "data_A_Vaga_PacBio/results/new_fragments_list.txt"
# matrixFile = "data_A_Vaga_PacBio/abs_fragments_contacts_weighted.txt"
# interactionFile = "data_A_Vaga_PacBio/interaction_matrix.pickle"
# outFile = "data_A_Vaga_PacBio/test_new_algo.gfa"

print('Loading the GFA file')
segments, names = load_gfa(gfaFile)
fragmentList = bf.read_fragment_list(fragmentsFile)

check_if_all_links_are_sorted(segments)

# Now computing the interaction matrix

interactionMatrix = bf.interactionMatrix(matrixFile, fragmentList, names, segments)

#interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))

#exporting it as to never have to do it again

print('Exporting interaction matrix')
file = open(interactionFile, 'wb')
pickle.dump(interactionMatrix, file)
    
file = open(interactionFile, 'rb')
interactionMatrix = pickle.load(file)



print("Solving ambiguities")

segments = solve_ambiguities(segments, interactionMatrix, names, 0.1, 0.2, 2)

print('Now exporting')

bf.export_to_GFA(segments, gfaFile = gfaFile, exportFile = outFile, merge_adjacent_contigs = False)
