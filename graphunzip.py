#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""

import input_output as io

# import analyse_HiC
from transform_gfa import gfa_to_fasta
from finish_untangling import merge_adjacent_contigs
from solve_with_long_reads import bridge_with_long_reads
from solve_with_HiC import solve_with_HiC
from determine_multiplicity import determine_multiplicity
#from segment import check_if_all_links_are_sorted

from scipy import sparse
import numpy as np
import argparse
import os.path
import sys
import pickle  # reading and writing files
import time


def parse_args_command() :
    
    parser = argparse.ArgumentParser()
    parser.add_argument("command", help="Either unzip, HiC-IM (to prepare Hi-C data) or linked-reads-IM (to prepare linked reads data)")
    
    return parser.parse_args(sys.argv[1:2])

def parse_args_unzip() :
    
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    groupInput = parser.add_argument_group("Input of GraphUnzip")
    groupOutput = parser.add_argument_group("Output of GraphUnzip")
    groupOther = parser.add_argument_group("Other options")
    
    groupInput.add_argument("-g", "--gfa", required=True, help="""GFA file to untangle""")
    groupInput.add_argument(
        "-i",
        "--HiCinteractions",
        required=False,
        default="Empty",
        help="""File containing the Hi-C interaction matrix from HiC-IM [default: None]""",
    )
    groupInput.add_argument(
        "-k",
        "--linkedReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the linked-reads interaction matrix from linked-reads-IM [default: None]""",
    )
    groupInput.add_argument(
        "-l", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format) or SPAligner (TSV format)"""
    )

    groupOutput.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output GFA [default: output.gfa]""",
    )
    groupOutput.add_argument(
        "-f",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
    )
    
    groupOther.add_argument(
        "-v",
        "--verbose",
        required = False,
        action="store_true",
    )
    # groupOther.add_argument(
    #     "-d",
    #     "--debug",
    #     required = False,
    #     default = '',
    #     help="""Activate the debug mode. Parameter: directory to put the logs and the intermediary GFAs.""",
    # )
    groupOther.add_argument(
        "--dont_merge",
        required=False,
        action="store_true",
        help="""If you don't want the output to have all possible contigs merged""",
    )
    
    groupOther.add_argument(
        "-u",
        "--unreliable_coverage",
        action="store_true",
        help="""Use this option if the coverage information of the graph is not reliable""",
    )
    
    return parser.parse_args(sys.argv[2:])

def parse_args_linked():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--gfa_graph", required=True,  help="""GFA file that will be untangled (required)""")
    
    parser.add_argument(
        "-p"
        "--linked_reads_IM", required=True, help = """Output file for the linked-read interaction matrix (required)""")
    
    parser.add_argument(
        "-b",
        "--barcoded_SAM", required=True, help = """SAM file of the barcoded reads aligned to the assembly. Barcodes must still be there (use option -C if aligning with BWA) (required)""")
    
    return parser.parse_args(sys.argv[2:])


def parse_args_HiC():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--gfa_graph", required=True,  help="""GFA file that will be untangled (required)""")
    
    parser.add_argument(
        "-m", "--matrix", required=True, help="""Sparse Hi-C contact map (required)"""
    )
    parser.add_argument(
        "-F", "--fragments", required=True, help="""Fragments list (required)"""
    )
    parser.add_argument(
        "--HiC_IM", required=False, default="Empty", help="""Output file for the Hi-C interaction matrix (required)"""
    )
    
    return parser.parse_args(sys.argv[2:])


def main():
    

    args_command = parse_args_command()
    command = args_command.command
    
    # if len(sys.argv) < 1 :
    #     sys.exit()
    
    t = time.time()
    
    if command == 'HiC-IM' :
        
        args = parse_args_HiC()
        
        matrixFile = args.matrix
        fragmentsFile = args.fragments
        
        outputIMH = args.HiC_IM
        
        gfaFile = args.gfa_graph
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
                    
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
            
    elif command == 'linked-reads-IM' :
        
        args = parse_args_linked()
        
        barcodedSAM = args.barcoded_SAM
        outputIMT = args.linked_reads_IM
        
        gfaFile = args.gfa_graph
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        
        if not os.path.exists(barcodedSAM):
            print('Error: could not find the SAM file.')
            sys.exit(1)
        
        tagInteractionMatrix = io.linkedReads_interactionMatrix(barcodedSAM, names)
        
        print("Exporting barcoded interaction matrix as ", outputIMT)
        with open(outputIMT, "wb") as o:
            pickle.dump(tagInteractionMatrix, o)

    elif command == 'unzip' :
        
        args = parse_args_unzip()
        
        gfaFile = args.gfa
        
        outFile = args.output
        fastaFile = args.fasta_output
        
        lrFile = args.longreads        
        
        interactionFileH = args.HiCinteractions
        interactionFileT = args.linkedReadsInteractions
        
        verbose = args.verbose
        # dbgDir = args.debug 
        merge = not args.dont_merge
        reliableCoverage = not args.unreliable_coverage
        
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        
        interactionMatrix = sparse.csr_matrix((len(segments), len(segments)))
        tagInteractionMatrix = sparse.csr_matrix((len(segments), len(segments)))
        useHiC = False
        uselr = False
        useTag = False
        
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

        print("Everything loaded, moving on to untangling the graph")
        
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
            segments = bridge_with_long_reads(segments, names, cn, lrFile, supported_links2, refHaploidy, multiplicities)
        
        #As a second step, use Hi-C and/or linked reads 
        if interactionMatrix.count_nonzero() > 0 :
            segments = solve_with_HiC(segments, interactionMatrix, names, confidentCoverage=reliableCoverage, verbose = verbose)
            print("Done untangling the graph")
        
        elif tagInteractionMatrix.count_nonzero() > 0 :
            segments = solve_with_HiC(segments, tagInteractionMatrix, names, confidentCoverage=reliableCoverage, verbose = verbose)

        elif not uselr :
            
            print("WARNING: all interaction matrices are empty, GraphUnzip does not do anything")
        
        print("Now exporting the result")
        merge_adjacent_contigs(segments)

            
        # now exporting the output
        io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge)
    
        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile)
    
        print("Finished in ", time.time() - t, " seconds")
        
    else :
        print("Unrecognized command ", command, "\". Use either unzip, HiC-IM (to prepare Hi-C data) or linked-reads-IM (to prepare linked reads data)")


if __name__ == "__main__":
    main()
