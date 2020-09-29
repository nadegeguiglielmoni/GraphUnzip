#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""

import input_output as io

# import analyse_HiC
from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities
from segment import check_if_all_links_are_sorted

from scipy import sparse
import numpy as np
import argparse
import os.path
import sys
import pickle  # reading and writing files
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
        "-fo",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
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
        default="Empty",
        help="""File with interactions [default: None]""",
    )
    parser.add_argument(
        "--merge",
        required=False,
        default="Empty",
        help="""If you want the output to have all possible contigs merged (y/n) [default: n]""",
    )
    return parser.parse_args()


def main():

    args = parse_args()
    gfaFile = args.gfa
    outFile = args.output
    fastaFile = args.fasta_output
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
    print("Loading the GFA file")
    segments, names = io.load_gfa(
        gfaFile
    )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...

    if fragmentsFile is not "Empty" and matrixFile is not "Empty":
        if os.path.exists(fragmentsFile):
            fragmentList = io.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix
            interactionMatrix = io.interactionMatrix(
                matrixFile, fragmentList, names, segments
            )

            if interactionFile is "Empty":
                interactionFile = "interactionMatrix.pickle"

            # exporting it as to never have to do it again

            print("Exporting interaction matrix")
            with open(interactionFile, "wb") as o:
                pickle.dump(interactionMatrix, o)

        else:
            print("Error: could not find fragments file {0}.".format(fragmentsFile))
            sys.exit(1)

    elif interactionFile is not "Empty":
        print("Loading the interaction matrix")
        interactionMatrix = io.load_interactionMatrix(interactionFile, segments, names)
    else:
        if not os.path.exists(interactionFile):
            print(
                "Error: you should provide either a processed interaction file, or the fragments list and the sparse contact map."
            )
            sys.exit(1)

    print("Everything loaded, moving on to solve_ambiguities")

    segments = solve_ambiguities(
        segments, interactionMatrix, names, stringenceReject, stringenceAccept, steps
    )

    # now exporting the output
    print("Now exporting")
    merge_adj = False
    if merge != "Empty" and merge != "n":
        merge_adj = True

    io.export_to_GFA(
        segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge_adj
    )

    if fastaFile != "None":
        gfa_to_fasta(outFile, fastaFile)

    print("Finished in ", time.time() - t, " seconds")


if __name__ == "__main__":
    main()

# gfaFile = "Arabidopsis/Arabidopsis_hybrid/simplified_graph.gfa"
# #gfaFile = "Arabidopsis/Arabidopsis_hybrid/small2.gfa"
# fragmentsFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/fragments_list.txt"
# matrixFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/abs_fragments_contacts_weighted.txt"
# interactionFile = "Arabidopsis/Arabidopsis_hybrid/unzip_out/interaction_matrix.pickle"
# outFile = "Arabidopsis/Arabidopsis_hybrid/unzip_out/unzipped.gfa"

# gfaFile = "data_A_Vaga_PacBio/Assembly.gfa"
# fragmentsFile = "data_A_Vaga_PacBio/mapping/new_fragments_list.txt"
# matrixFile = "data_A_Vaga_PacBio/mapping/abs_fragments_contacts_weighted.txt"
# interactionFile = "data_A_Vaga_PacBio/mapping/interaction_matrix.pickle"
# outFile = "data_A_Vaga_PacBio/unzipped.gfa"

# gfaFile = "Escherichia_Coli/1a1k/assemblyGraph_k63_noOverlaps.gfa"
# fragmentsFile = "Escherichia_Coli/1a1k/mapping/fragments_list.txt"
# matrixFile = "Escherichia_Coli/1a1k/mapping/abs_fragments_contacts_weighted.txt"
# interactionFile = "Escherichia_Coli/1a1k/mapping/interaction_matrix.pickle"
# outFile = "Escherichia_Coli/1a1k/unzipped.gfa"

# print('Loading the GFA file')
# segments, names = io.load_gfa(gfaFile)

# #check_if_all_links_are_sorted(segments)

# # Now computing the interaction matrix

# # fragmentList = bf.read_fragment_list(fragmentsFile)
# # interactionMatrix = bf.interactionMatrix(matrixFile, fragmentList, names, segments)
# # #interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))

# # # #exporting it as to never have to do it again

# # print('Exporting interaction matrix')
# # file = open(interactionFile, 'wb')
# # pickle.dump(interactionMatrix, file)

# #print(names)
# file = open(interactionFile, 'rb')
# interactionMatrix = io.load_interactionMatrix(file, segments, names)

# #print(names)
# #print(interactionMatrix[names['2276'], names['1504']])

# #print(interactionMatrix[names['111']])
# print(interactionMatrix[names['6003'], names['6002']])
# print(interactionMatrix[names['6004'], names['6002']])

# print('Next')

# #time.sleep(100)

# #print("Solving ambiguities")

# segments = solve_ambiguities(segments, interactionMatrix, names, stringenceReject = 0.1, stringenceAccept = 0.2, steps = 7)

# print('Now exporting')

# io.export_to_GFA(segments, gfaFile = gfaFile, exportFile = outFile, merge_adjacent_contigs = False)

# print('Done!')
