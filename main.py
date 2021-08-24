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
    
    parser.add_argument("command", help="Either unzip, HiC-IM, long-reads-IM or linked-reads-IM")
    
    groupUnzip = parser.add_argument_group("unzip options")
    groupHiC = parser.add_argument_group("HiC-IM options")
    grouplr = parser.add_argument_group("long-reads-IM options")
    grouplinked = parser.add_argument_group("linked-reads-IM options")
    
    parser.add_argument("-g", "--gfa", required = True, help="""GFA file to phase""")
    
    groupUnzip.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output GFA [default: output.gfa]""",
    )
    groupUnzip.add_argument(
        "-f",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
    )
    
    groupUnzip.add_argument(
        "-A",
        "--accepted",
        required=False,
        default=0.30,
        help="""Two links that are compared are deemed both true if
                        the weakest of the two, in term of Hi-C contacts, is
                        stronger than this parameter times the strength of the
                        strongest link [default: 0.30]""",
    )
    groupUnzip.add_argument(
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
    
    groupUnzip.add_argument(
        "-s",
        "--steps",
        required=False,
        default=10,
        help="""Number of cycles get rid of bad links - duplicate contigs. [default: 10]""",
    )

    groupHiC.add_argument(
        "-m", "--matrix", required=False, default="Empty", help="""Sparse Hi-C contact map"""
    )

    groupHiC.add_argument(
        "-F", "--fragments", required=False, default="Empty", help="""Fragments list"""
    )
    groupHiC.add_argument(
        "--HiC_IM", required=False, default="Empty", help="""Output file for the Hi-C interaction matrix (required)"""
    )
    
    groupUnzip.add_argument(
        "-i",
        "--HiCinteractions",
        required=False,
        default="Empty",
        help="""File containing the Hi-C interaction matrix from HiC-IM [default: None]""",
    )
    
    groupUnzip.add_argument(
        "-j",
        "--longReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the long-reads interaction matrix from long-reads-IM [default: None]""",
    )
    
    groupUnzip.add_argument(
        "-k",
        "--linkedReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the linked-reads interaction matrix from linked-reads-IM [default: None]""",
    )
    
    grouplr.add_argument(
        "-l", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format)"""
    )
    grouplr.add_argument(
        "--long_reads_IM", required = False, default="Empty", help="""Output file for the long-read interaction matrix (required)"""
    )
    
    groupUnzip.add_argument(
        "-e",
        "--exhaustive",
        action="store_true",
        help = "Removes all links not found in the GAF file (recommended if you have enough reads)",
    )
    
    grouplr.add_argument(
        "-M",
        "--minimum_match",
        required = False,
        default = 0,
        help = "Filters out alignments with a minimum match identity < minimum-match [default: 0]",
    )
    
    grouplr.add_argument(
        "-w",
        "--whole_match",
        action="store_true",
        help = "Filters out alignments that do not extend over the whole length of the read (recommended if you have enough reads)",
    )
    
    grouplinked.add_argument(
        "--linked_reads_IM", required=False, default = "Empty", help = """Output file for the linked-read interaction matrix (required)""")
    
    grouplinked.add_argument(
        "--barcoded_SAM", required=False, default = "Empty", help = """SAM file of the barcoded reads aligned to the assembly. Barcodes must still be there (use option -C if aligning with BWA) (required)""")
    
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

    
    groupUnzip.add_argument(
        "-v",
        "--verbose",
        required = False,
        action="store_true",
    )
    groupUnzip.add_argument(
        "-d",
        "--debug",
        required = False,
        default = '',
        help="""Activate the debug mode. Parameter: directory to put the logs and the intermediary GFAs.""",
    )
    groupUnzip.add_argument(
        "--merge",
        required=False,
        action="store_true",
        help="""If you want the output to have all possible contigs merged""",
    )
    return parser.parse_args()


def main():

    args = parse_args()
    
    command = args.command
    gfaFile = args.gfa
    
    outFile = args.output
    fastaFile = args.fasta_output
    
    matrixFile = args.matrix
    lrFile = args.longreads
    fragmentsFile = args.fragments
    barcodedSAM = args.barcoded_SAM
    
    outputIMH = args.HiC_IM
    outputIML = args.long_reads_IM
    outputIMT = args.linked_reads_IM
    
    interactionFileH = args.HiCinteractions
    interactionFileL = args.longReadsInteractions
    interactionFileT = args.linkedReadsInteractions
    
    stringenceReject = float(args.rejected)
    stringenceAccept = float(args.accepted)
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
    lrInteractionMatrix  = sparse.dok_matrix((len(segments), len(segments)))
    tagInteractionMatrix = sparse.dok_matrix((len(segments), len(segments)))
    lrLinks = []
    repeats = [] #this array will be uselful to flatten small loops based on long reads
    useHiC = False
    uselr = False
    useTag = False

    if command == 'HiC-IM' :
        
        if (fragmentsFile is not "Empty") and (matrixFile is not "Empty") and (outputIMH is not "Empty"):
            
            if os.path.exists(fragmentsFile) and os.path.exists(matrixFile):
                
                fragmentList = io.read_fragment_list(fragmentsFile)
    
                # Now computing the interaction matrix
    
                interactionMatrix = io.interactionMatrix(matrixFile, fragmentList, names, segments)
                useHiC = True
    
                # exporting it as to never have to do it again
    
                print("Exporting Hi-C interaction matrix as ", outputIMH)
                with open(outputIMH, "wb") as o:
                    pickle.dump(interactionMatrix, o)
    
            else:
                print("Error: could not find fragments file {0}.".format(fragmentsFile), " or matrix file {0}".format(matrixFile))
                sys.exit(1)    
        
        else :
            print("ERROR : options --fragments, --matrix and --HiC_IM are mandatory to use command HiC-IM")
     
    elif command == 'long-reads-IM' :
        
        if lrFile is not "Empty":
            
            if not os.path.exists(lrFile):
                print('Error: could not find the long-reads file.')
                sys.exit(1)
                
            lrInteractionMatrix, repeats, lrLinks = io.longReads_interactionsMatrix(lrFile, names, segments , similarity_threshold = mm, whole_mapping = wm)
            uselr = True
            
            if lrInteractionMatrix.count_nonzero() == 0 :
                print('WARNING: Tried loading the gaf file, but it seems that nothing could be read. If using mm and wm parameters, try lowering the minimum_match parameter towards 0 and/or the whole_match parameter to False. You can try to check the format of the gaf file and check if the names there correspond to the names of the GFA.')
                
            print("Exporting long-reads interaction matrix as ", outputIML)
            with open(outputIML, "wb") as o:
                pickle.dump(lrInteractionMatrix, o)
            with open(outputIML+'-r', "wb") as o:
                pickle.dump(repeats, o)
            with open(outputIML+'-l', "wb") as o:
                pickle.dump(lrLinks, o)
                
        else :
            
            print('ERROR: Providing a file with the long reads on the GFA (e.g. with GraphAligner) is mandatory to use the long-reads-IM, using option --longreads')
            
    elif command == 'linked-reads-IM' :
        
        if barcodedSAM is not "Empty" :
        
            if not os.path.exists(barcodedSAM):
                print('Error: could not find the SAM file.')
                sys.exit(1)
            
            tagInteractionMatrix = io.linkedReads_interactionMatrix(barcodedSAM, names)
            
            print("Exporting barcoded interaction matrix as ", outputIMT)
            with open(outputIMT, "wb") as o:
                pickle.dump(tagInteractionMatrix, o)
        
        else :
            print("ERROR: Providing the SAM of the barcoded reads aligned on the assembly is mandatory to use the linked-reads-IM command")
            
    elif command == 'unzip' :
        
        if interactionFileH is not "Empty":
            
            if not os.path.exists(interactionFileH) :
                print("ERROR: could not access ", interactionFileH)
                sys.exit(1)
            
            print("Loading the Hi-C interaction matrix")
            interactionMatrix = io.load_interactionMatrix(interactionFileH, segments, names, HiC = True)
            useHiC = True

        if interactionFileL is not "Empty" :
            
            if not os.path.exists(interactionFileL) :
                print("ERROR: could not access ", interactionFileL)
                sys.exit(1)
            
            print("Loading the long-reads interaction matrix")
            lrInteractionMatrix = io.load_interactionMatrix(interactionFileL, segments, names, HiC = False)
            uselr = True
            
            if not os.path.exists(interactionFileL+'-r') :
                print('WARNING: could not load ', interactionFileL+'-r', ', which should have been created along ', interactionFileL, '. Continuing since it is not absolutely necessary, but the program would work best with this file.')
                repeats = [0 for i in range(len(segments))]
            else :
                fi = open(interactionFileL+'-r', 'rb')
                repeats = pickle.load(fi)
                
            if not os.path.exists(interactionFileL+'-l') :
                print('WARNING: could not load ', interactionFileL+'-l', ', which should have been created along ', interactionFileL, '. This file is used for option -e, I turn it off and continue')
                exhaustive = False
            else :
                fi = open(interactionFileL+'-l', 'rb')
                lrLinks = pickle.load(fi)
                
        if interactionFileT is not "Empty":
            
            if not os.path.exists(interactionFileT) :
                print("ERROR: could not access ", interactionFileT)
                sys.exit(1)
            
            print("Loading the linked-reads interaction matrix")
            tagInteractionMatrix = io.load_interactionMatrix(interactionFileT, segments, names, HiC = False)
            useTag = True
            
        if not( useHiC or uselr or useTag) :
            
            print("ERROR: You should provide to unzip interaction matrices, using either --HiCinteractions (-i), --longReadsInteractions (-j) or --linkedReadsInteractions (-k). If you do not have them, you can create them using the HiC-IM, long-reads-IM or linked-reads-IM commands")


        print("Everything loaded, moving on to unzipping")
        cn = {}
        
        if interactionMatrix.count_nonzero() > 0 or lrInteractionMatrix.count_nonzero() >0 or tagInteractionMatrix.count_nonzero() > 0 or exhaustive:
                    
            segments, cn = solve_ambiguities(
                segments, interactionMatrix, lrInteractionMatrix, tagInteractionMatrix, names, stringenceReject, stringenceAccept, steps, repeats = repeats, copiesNumber = cn, debugDir = dbgDir, lr_links = lrLinks, check_links = exhaustive, verbose = verbose,
            )
        
        else :
            
            print("WARNING: all interaction matrices are empty, GraphUnzip does not do anything")
            
        if lrInteractionMatrix.count_nonzero() == 0 and uselr:
            print("WARNING: the long reads interaction matrix between contigs is empty. This could be due to having filtered out all information from long reads. If you used --exhaustive I remove all edges, I do nothing elsewhise.")
    
        # now exporting the output
        print("Now exporting")
    
        io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge)
    
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
