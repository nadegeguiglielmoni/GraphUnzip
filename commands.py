#!/usr/bin/python
# -*- coding: utf-8 -*-
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example

from docopt import docopt

import input_output as io

from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities

# from segment import check_if_all_links_are_sorted

from scipy import sparse
import numpy as np
import os.path
import sys
import pickle  # reading and writing files
import time


class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """Initialize the commands."""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError


class GFAtoFasta(AbstractCommand):
    """ GFA to fasta command
    Convert GFA file to fasta file.

    usage:
        gfa2fasta [--outfile=genome.fasta] <genome.gfa>

    arguments:
        genome.gfa            Draft assembly in GFA format.
        
    options:
        -o, --outfile=FILE    Fasta output [default: genome.fasta].
    """

    def execute(self):
        print("Running gfa2fasta module.")
        gfa_to_fasta(self.args["<mapping.bam>"], self.args["--outfile"])


class ProcessHiCData(AbstractCommand):
    """ Processing Hi-C data command
    Process Hi-C sparse matrix for GraphUnzip.

    usage:
        processHiCData [--outfile=hicMatrix.pickle] --gfa=FILE --frags=FILE 
                       <sparse_matrix.txt>

    arguments:
        sparse_matrix.txt     Sparse Hi-C matrix.
        
    options:
        -o, --outfile=FILE    Processed Hi-C data output [default: hicMatrix.pickle].
        -g, --gfa=FILE        Draft assembly in GFA format.
        -F, --frags=FILE      Fragments list.
    """

    def execute(self):
        print("Running Hi-C data processing module.")

        print("Loading GFA file.")
        segments, names = io.load_gfa(file=self.args["--gfa"])

        print("Loading fragment list.")
        fragmentList = io.read_fragment_list(file=self.args["--frags"], header=True)

        hicInteractionMatrix = io.interactionMatrix(
            hicFile=self.args["<sparse_matrix.txt>"],
            fragmentList=fragmentList,
            names=names,
            segments=segments,
            header=True,
        )

        print("Exporting interaction matrix as ", self.args["--outfile"])
        with open(self.args["--outfile"], "wb") as o:
            pickle.dump(hicInteractionMatrix, o)


class ProcessLRData(AbstractCommand):
    """ Processing long-read data command
    Process mapped long-read data for GraphUnzip.

    usage:
        processLRData [--outfile=LRMatrix.pickle] [--outlinks=links.pickle] 
                      [--min-match=0] [--whole-match] --gfa=FILE <mapping.gaf>

    arguments:
        mapping.gaf             Long reads mapped to the GFA.
        
    options:
        -o, --outfile=FILE      Processed long-read data output [default: LRMatrix.pickle].
        -l, --outlinks=FILE     Links output [default: links.pickle].
        -m, --min-match=INT     Minimum match threshold [default: 0].
        -w, --whole-match       Keep only alignments of full reads.
        -g, --gfa=FILE          Draft assembly in GFA format.
    """

    def execute(self):
        print("Running long-read data processing module.")

        print("Loading GFA file.")
        segments, names = io.load_gfa(file=self.args["--gfa"])

        lrInteractionMatrix, lrLinks = io.longReads_interactionsMatrix(
            gafFile=self.args["<mapping.gaf>"],
            names=names,
            segments=segments,
            similarity_threshold=self.args["--min-match"],
            whole_mapping=self.args["--whole-match"],
        )

        print("Exporting interaction matrix as ", self.args["--outfile"])
        with open(self.args["--outfile"], "wb") as o:
            pickle.dump(lrInteractionMatrix, o)

        print("Exporting links as ", self.args["--outlinks"])
        with open(self.args["--outlinks"], "wb") as o:
            pickle.dump(lrLinks, o)


class Unzip(AbstractCommand):
    """ Unzip command
    Unzip assembly graph.

    usage:
        unzip [--hic=hicMatrix.pickle] [--lr=LRMatrix.pickle] [--lr-links=links.pickle] 
              [--min-match=0] [--whole-match] --gfa=FILE <mapping.gaf>

    arguments:
        mapping.gaf            Long reads mapped to the GFA.
        
    options:
        -h, --hic=FILE         Processed Hi-C data [default: hicMatrix.pickle].
        -l, --lr=FILE          Processed long-read data [default: LRMatrix.pickle].
        -L, --lr-links=FILE    Links output [default: links.pickle].
        -g, --gfa=FILE         Draft assembly in GFA format.
        -A, --accept           Threshold to accept links [default: 0.30].
        -R, --reject           Threshold to reject links [default: 0.15].
        -S, --steps            Number of cycles to get rid of bad links [default: 10].
        -e, --exhaustive       Remove all links not found in the GAF file.
    """

    def execute(self):
        print("Running unzip module.")

        print("Loading GFA file.")
        segments, names = io.load_gfa(file=self.args["--gfa"])

        print("Loading Hi-C interaction matrix.")
        interactionMatrix = io.load_interactionMatrix(interactionFile, segments, names)

        if interactionMatrix.count_nonzero() > 0 or exhaustive:
            segments, cn = solve_ambiguities(
                segments,
                interactionMatrix,
                names,
                self.args["--reject"],
                self.args["--accept"],
                steps,
                SEGMENT_REPEAT=normalizationFactor * 10,
                copiesNumber=cn,
                debugDir=dbgDir,
                lr_links=lrLinks,
                check_links=exhaustive,
            )
