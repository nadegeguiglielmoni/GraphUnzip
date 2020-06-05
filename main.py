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

from copy import deepcopy

import argparse
import os.path
import sys

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
        default=0.45,
        help="""Threshold to accept links. [default: 0.40]""",
    )
    parser.add_argument(
        "-R",
        "--rejected",
        required=False,
        default=0.20,
        help="""Threshold to reject links. [default: 0.20]""",
    )
    parser.add_argument(
        "-s",
        "--steps",
        required=False,
        default=10,
        help="""Number of cycles get rid of bad links - duplicate contigs. [default: 10]""",
    )
    # parser.add_argument(
    #     "-f",
    #     "--fasta",
    #     required=False,
    #     default="Empty",
    #     help="""Segments from the GFA in fasta format""",
    # )
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

    if not os.path.exists(gfaFile):
        print("Error: could not find GFA file {0}.".format(gfaFile))
        sys.exit(1)

    #Creating/loading the interaction matrix
    if fragmentsFile is not "Empty" and matrixFile is not "Empty":
        if os.path.exists(fragmentsFile):
            fragmentList = bf.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix
            interactionMatrix = bf.interactionMatrix(matrixFile, fragmentList)
            print('Interaction matrix built')

            # exporting it as to never have to do it again
            # if not os.path.exists(interactionFile):
            #     bf.export_to_csv(interactionMatrix, interactionFile)
            # else:
            #     print(
            #         "Error: {0} already exists, please remove it.".format(
            #             interactionFile
            #         )
            #     )
            #     sys.exit(1)
        else:
            print("Error: could not find fragments file {0}.".format(fragmentsFile))
            sys.exit(1)
    else:
        if not os.path.exists(interactionFile):
            print("Error: you should provide either a processed interaction file, or the fragments list and the sparse contact map.")
            sys.exit(1)
            
    # Loading the data
    originalLinks, CIGARlinks, names, lengths = load_gfa(gfaFile)

    #interactionMatrix = bf.import_from_csv(interactionFile)
    
    print('Everything loaded, beginning to solve_ambiguities')
    links, listOfSuperContigs, copiesNumber = solve_ambiguities(deepcopy(originalLinks), names, interactionMatrix, lengths, lambda x:1, stringenceReject, stringenceAccept, steps)

    # now exporting the output
    bf.export_to_GFA(links, listOfSuperContigs, copiesNumber, originalLinks, CIGARlinks, names, gfaFile, exportFile=outFile)


if __name__ == "__main__":
    main()


# originalLinks, CIGARlinks, names, lengths = load_gfa('data_A_Vaga_PacBio/Assembly.gfa')
# interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')

# print('Loaded')
 
# links, listOfSuperContigs, cn = solve_ambiguities(deepcopy(originalLinks), names, interactionMatrix, lengths, lambda x:1, 0.2, 0.45 ,15) #rejectedThreshold<AcceptedThreshold
# bf.export_to_GFA(links, listOfSuperContigs, cn, originalLinks, originalLinksCIGAR = CIGARlinks, names = names, gfaFile = 'data_A_Vaga_PacBio/Assembly.gfa', exportFile = 'results/A_Vaga_PacBio/A_Vaga_finished2.gfa')