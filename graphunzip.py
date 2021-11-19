#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""

import input_output as io

# import analyse_HiC
from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities
from solve_with_long_reads import bridge_with_long_reads
from determine_multiplicity import determine_multiplicity
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
        "-k",
        "--linkedReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the linked-reads interaction matrix from linked-reads-IM [default: None]""",
    )
    
    groupUnzip.add_argument(
        "-l", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format), if you have them"""
    )
    
    # groupUnzip.add_argument(
    #     "-e",
    #     "--exhaustive",
    #     action="store_true",
    #     help = "Removes all links not found in the GAF file (recommended if you have enough reads)",
    # )
    
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
    #     "-Rlr",git c
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
        "--dont_merge",
        required=False,
        action="store_true",
        help="""If you don't want the output to have all possible contigs merged""",
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
    outputIMT = args.linked_reads_IM
    
    interactionFileH = args.HiCinteractions
    interactionFileT = args.linkedReadsInteractions
    
    stringenceReject = float(args.rejected)
    stringenceAccept = float(args.accepted)
    steps = int(args.steps)
    
    verbose = args.verbose
    
    dbgDir = args.debug
    
    merge = not args.dont_merge

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
    tagInteractionMatrix = sparse.dok_matrix((len(segments), len(segments)))
    useHiC = False
    uselr = False
    useTag = False

    if command == 'HiC-IM' :
        
        if (fragmentsFile != "Empty") and (matrixFile != "Empty") and (outputIMH != "Empty"):
            
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
            
    elif command == 'linked-reads-IM' :
        
        if barcodedSAM != "Empty" :
        
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
        
        if interactionFileH != "Empty":
            
            if not os.path.exists(interactionFileH) :
                print("ERROR: could not access ", interactionFileH)
                sys.exit(1)
            
            print("Loading the Hi-C interaction matrix")
            interactionMatrix = io.load_interactionMatrix(interactionFileH, segments, names, HiC = True)
            useHiC = True

        if lrFile != "Empty" :
            
            if not os.path.exists(lrFile) :
                print("ERROR: could not access ", lrFile)
                sys.exit(1)   
            uselr = True
            
                
        if interactionFileT != "Empty":
            
            if not os.path.exists(interactionFileT) :
                print("ERROR: could not access ", interactionFileT)
                sys.exit(1)
            
            print("Loading the linked-reads interaction matrix")
            tagInteractionMatrix = io.load_interactionMatrix(interactionFileT, segments, names, HiC = False)
            useTag = True
            
        if not( useHiC or uselr or useTag) :
            
            print("ERROR: You should provide to unzip long reads mapped in GAF format and/or interaction matrices, using either --HiCinteractions (-i) or --linkedReadsInteractions (-k). If you do not have them, you can create them using the HiC-IM, long-reads-IM or linked-reads-IM commands")
            sys.exit()

        print("Everything loaded, moving on to unzipping")
        
        #creating copiesnuber (cn), a dictionnary inventoring how many times 
        cn = {}
        for segment in segments :
            for name in segment.names :
                cn[name] = 1
        
        ##Moving to the actual unzipping of the graph
        
        supported_links2 = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links considering the topography of the graph
        refHaploidy, multiplicities = determine_multiplicity(segments, names, supported_links2) #multiplicities can be seen as a mininimum multiplicity of each contig regarding the topology of the graph
        
        #As a first step, use only the long reads, if available
        if uselr :
            bridge_with_long_reads(segments, names, cn, lrFile, supported_links2, refHaploidy, multiplicities)
        
        #As a second step, use Hi-C and/or linked reads 
        if interactionMatrix.count_nonzero() > 0 or tagInteractionMatrix.count_nonzero() > 0 :
                    
            segments, cn = solve_ambiguities(segments, interactionMatrix, tagInteractionMatrix, multiplicities, names, stringenceReject, stringenceAccept, steps, copiesNumber = cn, debugDir = dbgDir, verbose = verbose)
        
        elif not uselr :
            
            print("WARNING: all interaction matrices are empty, GraphUnzip does not do anything")
            
            
        # now exporting the output
        print("Now exporting")
    
        io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge)
    
        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile)
    
        print("Finished in ", time.time() - t, " seconds")


if __name__ == "__main__":
    main()
