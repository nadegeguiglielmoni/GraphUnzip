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
import argparse
import os.path
import sys


def parse_args():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gfa", required=True, help="""GFA file""")
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
        help="""Number of steps to get rid of bad links. [default: 10]""",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=False,
        default="Empty",
        help="""Segments from the GFA in fasta format""",
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
    return parser.parse_args()


def main():

    args = parse_args()
    gfaFile = args.gfa
    fastaFile = args.fasta
    matrixFile = args.matrix
    fragmentsFile = args.fragments
    interactionFile = args.interaction
    stringenceReject = float(args.rejected)
    stringenceAccept = float(args.accepted)
    steps = int(args.step)

    if not os.path.exists(gfaFile):
        print("Error: could not find GFA file {0}.".format(gfaFile))
        sys.exit(1)

    # Loading the data
    links, names = load_gfa(gfaFile)

    if fragmentsFile is not "Empty" and matrixFile is not "Empty":
        if os.path.exists(fragmentsFile):
            fragmentList = bf.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix
            interactionMatrix = bf.interactionMatrix(matrixFile, fragmentList)

            # exporting it as to never have to do it again
            if not os.path.exists(interactionFile):
                bf.export_to_csv(interactionMatrix, interactionFile)
            else:
                print(
                    "Error: {0} already exists, please remove it.".format(
                        interactionFile
                    )
                )
                sys.exit(1)
        else:
            print("Error: could not find fragments file {0}.".format(fragmentsFile))
            sys.exit(1)
    else:
        if not os.path.exists(interactionFile):
            print(
                "Error: you should provide either a processed interaction file, or the fragments list and the sparse contact map."
            )
            sys.exit(1)

    if fastaFile is "Empty":
        gfa_to_fasta(gfaFile, "gfa_to_fasta.fasta")

    interactionMatrix = bf.import_from_csv("interactionMatrix.csv")

    links, listOfSuperContigs, copiesNumber = solve_ambiguities(
        links, interactionMatrix, stringenceReject, stringenceAccept, steps, names
    )

    # now exporting the output
    bf.export_to_GFA(
        links, listOfSuperContigs, copiesNumber, names, fastaFile, exportFile=outFile
    )


if __name__ == "__main__":
    main()
