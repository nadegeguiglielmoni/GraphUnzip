#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""

import input_output as io

# import analyse_HiC
from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities
#from segment import check_if_all_links_are_sorted

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
        default=0.30,
        help="""Two links that are compared are deemed both true if
                        the weakest of the two, in term of Hi-C contacts, is
                        stronger than this parameter times the strength of the
                        strongest link [default: 0.30]""",
    )
    parser.add_argument(
        "-R",
        "--rejected",
        required=False,
        default=0.15,
        help="""When two links are compared, the weakest of the two,
                        in term of Hi-C contacts, is considered false and
                        deleted if it is weaker than this parameter times the
                        strength of the strongest links (always smaller than
                        --accepted)[default: 0.15]""",
    )
    
    parser.add_argument(
        "-s",
        "--steps",
        required=False,
        default=10,
        help="""Number of cycles get rid of bad links - duplicate contigs. [default: 10]""",
    )

    parser.add_argument(
        "-m", "--matrix", required=False, default="Empty", help="""Sparse Hi-C contact map"""
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
        "-lr", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format)"""
    )
    
    parser.add_argument(
        "-e",
        "--exhaustive",
        action="store_true",
        help = "Removes all links not found in the GAF file (recommended if you have enough reads)",
    )
    
    parser.add_argument(
        "-mm",
        "--minimum_match",
        required = False,
        default = 0,
        help = "Filters out alignments with a minimum match identity < minimum-match [default: 0]",
    )
    
    parser.add_argument(
        "-wm",
        "--whole_match",
        action="store_true",
        help = "Filters out alignments that do not extend over the whole length of the read (recommended if you have enough reads)",
    )
    
    # parser.add_argument(
    #     "-Alr",
    #     "--accepted-lr",
    #     required=False,
    #     default=0.30,
    #     help="""Threshold to accept long-reads links. [default: 0.30]""",
    # )
    # parser.add_argument(
    #     "-Rlr",
    #     "--rejected-lr",
    #     required=False,
    #     default=0.15,
    #     help="""Threshold to reject long-reads links. [default: 0.15]""",
    # )
    # parser.add_argument(
    #     "-c",
    #     "--combine-matrices",
    #     required=False,
    #     default=0,
    #     help="""Define how Hi-C and long reads interaction matrices are combined. Values are 0 (add the two matrices), 1 (run first long-reads), 2 (run first Hi-C). [default: 0]""",
    # )
    
    parser.add_argument(
        "-v",
        "--verbose",
        required = False,
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--debug",
        required = False,
        default = '',
        help="""Activate the debug mode. Parameter: directory to put the logs and the intermediary GFAs.""",
    )
    parser.add_argument(
        "--merge",
        required=False,
        action="store_true",
        help="""If you want the output to have all possible contigs merged""",
    )
    return parser.parse_args()


def main():

    args = parse_args()
    
    gfaFile = args.gfa
    outFile = args.output
    fastaFile = args.fasta_output
    matrixFile = args.matrix
    lrFile = args.longreads
    fragmentsFile = args.fragments
    interactionFile = args.interactions
    stringenceReject = float(args.rejected)
    stringenceAccept = float(args.accepted)
    # stringenceRejectLR = float(args.rejected-lr)
    # stringenceAcceptLR = float(args.accepted-lr)
    steps = int(args.steps)
    exhaustive = bool(args.exhaustive)
    
    mm = float(args.minimum_match)
    wm = bool(args.whole_match)
    
    verbose = args.verbose
    
    dbgDir = args.debug
    
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

    interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))
    useHiC = False
    uselr = False


    if fragmentsFile is not "Empty" and matrixFile is not "Empty":
        if os.path.exists(fragmentsFile):
            fragmentList = io.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix

            interactionMatrix = io.interactionMatrix(
                matrixFile, fragmentList, names, segments
            )
            useHiC = True

            if interactionFile is "Empty":
                interactionFile = "interactionMatrix.pickle"

            # exporting it as to never have to do it again

            print("Exporting interaction matrix as ", interactionFile)
            with open(interactionFile, "wb") as o:
                pickle.dump(interactionMatrix, o)

        else:
            print("Error: could not find fragments file {0}.".format(fragmentsFile))
            sys.exit(1)    
        
    elif interactionFile is not "Empty":
        print("Loading the interaction matrix")
        interactionMatrix = io.load_interactionMatrix(interactionFile, segments, names)
        useHiC = True
        
    elif not os.path.exists(interactionFile) and not os.path.exists(lrFile):
        print(
            "Error: could not find the file(s) to build/load the interaction matrix. You should provide either a processed interaction file in pickle format or a fragment list and a hic contact sparse matrix in hicstuff format. You can check you spelled everything correctly."
        )
        sys.exit(1)
    
    lrLinks = []
    repeats = []
    lrInteractionMatrix  = sparse.dok_matrix((len(segments), len(segments)))
    normalizationFactor = 1
    if lrFile is not "Empty":
        
        if not os.path.exists(lrFile):
            print('Error: could not find the long-reads file.')
            sys.exit(1)
            
        lrInteractionMatrix, repeats, lrLinks = io.longReads_interactionsMatrix(lrFile, names, segments , similarity_threshold = mm, whole_mapping = wm)
        
        if lrInteractionMatrix.count_nonzero() == 0 :
             print('WARNING: Tried loading the gaf file, but it seems that nothing could be read. If using mm and wm parameters, try lowering the minimum_match parameter towards 0 and/or the whole_match parameter to False. You can try to check the format of the gaf file, check if the names there correspond to the names of the GFA. ')

        uselr = True
        # lrSum = np.sum([np.sum(i) for i in lrInteractionMatrix])
        # hicSum = np.sum([np.sum(i) for i in interactionMatrix])
        # normalizationFactor = hicSum/lrSum * 0.1 + 1 #you can use a factor bigger than 1 to give more importance to long reads, smaller than 1 to give more importance to Hi-C
        # 
        # for i in lrInteractionMatrix.keys() :
        #         
        #         lrInteractionMatrix[i] *= normalizationFactor
        #         print('a, ', i)
        
        #interactionMatrix += lrInteractionMatrix
    else :
        exhaustive = False # checking all links makes no sense if there are no links
    

    print("Everything loaded, moving on to solve_ambiguities")
    cn = {}
    
    if interactionMatrix.count_nonzero() > 0 or lrInteractionMatrix.count_nonzero() >0 or exhaustive:
        segments, cn = solve_ambiguities(
            segments, interactionMatrix, lrInteractionMatrix, names, stringenceReject, stringenceAccept, steps, repeats = repeats, copiesNumber = cn, debugDir = dbgDir, lr_links = lrLinks, check_links = exhaustive, verbose = verbose,
        )
    if lrInteractionMatrix.count_nonzero() == 0 and uselr:
        print("WARNING: the long reads interaction matrix between contigs is empty. This could be due to having filtered out all information from long reads. If you used --exhaustive I remove all edges, I do nothing elsewhise.")

    # now exporting the output
    print("Now exporting")

    io.export_to_GFA(
        segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge
    )

    if fastaFile != "None":
        io.export_to_fasta(segments, gfaFile, fastaFile)

    print("Finished in ", time.time() - t, " seconds")


if __name__ == "__main__":
    main()


# gfaFile = "A_vaga_article/Nanopore_Ratatosk/avaga.flye_keep-haplotypes_hifi.ont_ratatosk_all.gfa"
# interactionFile = "A_vaga_article/Nanopore_Ratatosk/ont_ratatosk_hicMatrix.pickle"

# gafFile = '../users/Cyril_Matthey-Doret/ont_to_gfa.gaf'
# gfaFile = '../users/Cyril_Matthey-Doret/assembly_graph.gfa'
# outFile =  '../users/Cyril_Matthey-Doret/output.gfa'
# 
# 
# print('Loading the GFA file')
# segments, names = io.load_gfa(gfaFile)

#check_if_all_links_are_sorted(segments)

#Now computing the interaction matrix

# lrinteractionMatrix, repeats, allLinks = io.longReads_interactionsMatrix(gafFile, names, segments, 0.9, True)
# interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))

# fragmentList = io.read_fragment_list(fragmentsFile)
# interactionMatrix = io.interactionMatrix(matrixFile, fragmentList, names, segments)
#interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))

# #exporting it as to never have to do it again

# print('Exporting interaction matrix')
# file = open(interactionFile, 'wb')
# pickle.dump(interactionMatrix, file)

#print(names)
# hicinteractionMatrix = io.load_interactionMatrix(interactionFile, segments, names)

# interactionMatrix = lrinteractionMatrix #+ hicinteractionMatrix

# #print(allLinks)

#print("Solving ambiguities")

# segments, cn = solve_ambiguities(segments, interactionMatrix, lrinteractionMatrix, names, stringenceReject = 0.15, stringenceAccept = 0.3, steps = 5, lr_links = allLinks, useNeighborOfNeighbor = False, debugDir = '../users/Cyril_Matthey-Doret/debug', check_links = True)
# 
# io.export_to_GFA(segments, gfaFile = gfaFile, exportFile = outFile+'.temp', merge_adjacent_contigs = False)

# # segments = solve_ambiguities(segments, lrinteractionMatrix, names, stringenceReject = 0.2, stringenceAccept = 0.4, steps = 7, SEGMENT_REPEAT = 10)


# print('Done!')
