#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 6 07:42:14 2020

"""

__version__="1.0.0"
__author__="Roland Faure"

import input_output as io

# import analyse_HiC
from transform_gfa import gfa_to_fasta
from finish_untangling import merge_adjacent_contigs
from solve_with_long_reads import bridge_with_long_reads
#from solve_with_long_reads2 import bridge_with_long_reads2
from solve_with_HiC import solve_with_HiC
from determine_multiplicity import determine_multiplicity
from clean_graph import clean_graph
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
    
    parser.add_argument("command", help="Either unzip (untangle the GFA file), purge (retain only haploid contigs), extract (extract haploid assembly with a close reference genome), HiC-IM (to prepare Hi-C data) or linked-reads-IM (to prepare linked reads data)")
    
    return parser.parse_args(sys.argv[1:2])

def parse_args_purge():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--gfa", required=True,  help="""GFA file to be purged (required)""")

    parser.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output the purged GFA assemby [default: output.gfa]""",
    )
    parser.add_argument(
        "-f",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
    )
    parser.add_argument(
        "--dont_merge",
        required=False,
        action="store_true",
        help="""If you don't want the output to have all possible contigs merged""",
    )
    
    return parser.parse_args(sys.argv[2:])


def parse_args_extract():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--gfa", required=True,  help="""GFA file from which the assembly will be extracted (required)""")
    
    parser.add_argument(
        "-l",
        "--genome", required=True, help="""Genome mapped to the GFA with GraphAligner (GAF format) or SPAligner (TSV format) (required)""")
    
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output the GFA assemby of the genome [default: output.gfa]""",
    )
    parser.add_argument(
        "-f",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
    )
    
    parser.add_argument(
        "-r",
        "--dont_rename",
        action="store_true",
        help="""Use if you don't want to name the resulting supercontigs with short names but want to keep the names of the original contigs""")
    
    parser.add_argument(
        "--dont_merge",
        required=False,
        action="store_true",
        help="""If you don't want the output to have all possible contigs merged""",
    )
    
    return parser.parse_args(sys.argv[2:])

def parse_args_unzip() :
    
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    groupInput = parser.add_argument_group("Input of GraphUnzip")
    groupOutput = parser.add_argument_group("Output of GraphUnzip")
    groupBehavior = parser.add_argument_group("Behavior of GraphUnzip")
    groupOther = parser.add_argument_group("Other options")
    
    
    groupInput.add_argument("-g", "--gfa", required = True, help="""GFA file to phase""")
 
    groupInput.add_argument(
        "-i",
        "--HiCinteractions",
        required=False,
        default="Empty",
        help="""File containing the Hi-C interaction matrix from HiC-IM [optional]""",
    )
    groupInput.add_argument(
        "-k",
        "--linkedReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the linked-reads interaction matrix from linked-reads-IM [optional]""",
    )
    groupInput.add_argument(
        "-l", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format) or SPAligner (TSV format) [optional]"""
    )

    groupInput.add_argument(
        "-s", "--genomeSize", required = False, default="Empty", help="""Full genome size, counting all haplotypes - e.g. 100m or 3g [optional but recommended]"""
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
    groupOutput.add_argument(
        "-b",
        "--bam_file",
        required=False,
        default="None",
        help="""bam file of the Hi-C reads aligned on assembly. GraphUnzip will output bam_file.new.bam corresponding to the new bam file, ready to be used for scaffolding [optional]""",
    )
    
    groupOther.add_argument(
        "-v",
        "--verbose",
        required = False,
        action="store_true",
    )
    groupOther.add_argument(
        "-r",
        "--dont_rename",
        action="store_true",
        help="""Use if you don't want to name the resulting supercontigs with short names but want to keep the names of the original contigs""",
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
    
    groupBehavior.add_argument(
        "-H",
        "--haploid",
        action="store_true",
        help="""Use this option if you wish to obtain a collapsed assembly of a multiploid genome.""",
    )
    groupBehavior.add_argument(
        "-c",
        "--conservative",
        action="store_true",
        help="""(Hi-C only) Output very robust contigs. Use this option if the coverage information of the graph is not reliable""",
    )
    
    groupBehavior.add_argument(
        "-B",
        "--bold",
        action="store_true",
        help="""(Hi-C only)[default] Proposes the best untangling it can get (can be misled by approximate coverage information). Use this option if the contig coverage information of the graph can be trusted""",
    )
    groupBehavior.add_argument(
        "-n",
        "--noisy",
        action="store_true",
        help="""(Hi-C only) Use this option if you expect that the assembly may contain artefactual contigs, e.g. when you use the .p_utg.gfa of hifiasm""",
    )
    groupBehavior.add_argument(
        "-e",
        "--exhaustive",
        action="store_true",
        help="""(long reads only) All links not found in the .gaf will be removed""",
    )

    
    return parser.parse_args(sys.argv[2:])

def parse_args_linked():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--gfa_graph", required=True,  help="""GFA file that will be untangled (required)""")
    
    parser.add_argument(
        "-p",
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
        "-b", "--bam", default="Empty", required=False, help="""Bam file of Hi-C reads aligned on assembly and sorted by name (if using bam format)"""
    )
    parser.add_argument(
        "-m", "--matrix", default="Empty", required=False, help="""Sparse Hi-C contact map (if using instaGRAAL format)"""
    )
    parser.add_argument(
        "-F", "--fragments", default="Empty", required=False, help="""Fragments list (if using instaGRAAL format)"""
    )
    parser.add_argument(
        "-i","--HiC_IM", required=True, help="""Output file for the Hi-C interaction matrix (required)"""
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
        bamfile = args.bam
        
        if bamfile == "Empty" and (matrixFile == "Empty" or fragmentsFile == "Empty"):
            print("ERROR: you must provide as input either (a bam file) or (an abs_fragments_weighted.txt and a fragment_list.txt files from hicstuff)")
            sys.exit()
        
        outputIMH = args.HiC_IM
        
        gfaFile = args.gfa_graph
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0 :
            print("ERROR: could not read the GFA")
            sys.exit()
        
        if bamfile == "Empty":
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
        else:
            if os.path.exists(bamfile):
                interactionMatrix = io.read_bam(bamfile, names, segments)
                useHiC = True
                
                print("Exporting Hi-C interaction matrix as ", outputIMH)
                with open(outputIMH, "wb") as o:
                    pickle.dump(interactionMatrix, o)  
            else:
                print("Error: could not find bam file {0}.".format(bamfile))
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
        if len(segments) == 0 :
            print("ERROR: could not read the GFA")
            sys.exit()
        
        if not os.path.exists(barcodedSAM):
            print('Error: could not find the SAM file.')
            sys.exit(1)
        
        tagInteractionMatrix = io.linkedReads_interactionMatrix(barcodedSAM, names)
        
        print("Exporting barcoded interaction matrix as ", outputIMT)
        with open(outputIMT, "wb") as o:
            pickle.dump(tagInteractionMatrix, o)
            
    elif command == 'extract' :
        
        args = parse_args_extract()

        gfaFile = args.gfa
        outFile = args.output
        fastaFile = args.fasta_output
        
        lrFile = args.genome
        
        rename = not args.dont_rename
        merge = not args.dont_merge
        
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0 :
            print("ERROR: could not read the GFA")
            sys.exit()
            
        #creating copiesnuber (cn), a dictionnary inventoring how many times each contig appears
        cn = {}
        for segment in segments :
            for name in segment.names :
                cn[name] = 1
        
        print("================\n\nEverything loaded, moving on to untangling the graph\n\n================")
        
        supported_links2 = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links considering the topography of the graph
        refHaploidy, multiplicities = determine_multiplicity(segments, supported_links2, reliable_coverage=False) #multiplicities can be seen as a mininimum multiplicity of each contig regarding the topology of the graph

        segments = bridge_with_long_reads(segments, names, cn, lrFile, supported_links2, multiplicities, exhaustive=True, extract=True)
        print("Merging contigs that can be merged...")
        merge_adjacent_contigs(segments)
        print("\n*Done extracting the genome*\n")
        
        # now exporting the output  
        print("Now exporting the result")
        io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge, rename_contigs=rename)
    
        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile, rename_contigs=rename)
    
        print("Finished in ", time.time() - t, " seconds")

    elif command == 'unzip' :
        
        args = parse_args_unzip()
        
        gfaFile = args.gfa
        
        outFile = args.output
        fastaFile = args.fasta_output
        bamFile = args.bam_file
        
        lrFile = args.longreads        
        
        interactionFileH = args.HiCinteractions
        interactionFileT = args.linkedReadsInteractions

        
        verbose = args.verbose
        rename = not args.dont_rename
        # dbgDir = args.debug 
        haploid = args.haploid
        merge = not args.dont_merge
        reliableCoverage = not args.conservative
        exhaustive = args.exhaustive
        noisy = args.noisy

        genomeSize = 0 
        #convert genome size to integer
        if args.genomeSize != "Empty":
            if args.genomeSize[-1] == "g":
                genomeSize = int(args.genomeSize[:-1])*1000000000
            elif args.genomeSize[-1] == "m":
                genomeSize = int(args.genomeSize[:-1])*1000000
            elif args.genomeSize[-1] == "k":
                genomeSize = int(args.genomeSize[:-1])*1000
            else:
                genomeSize = int(args.genomeSize)
                
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0 :
            print("ERROR: could not read the GFA")
            sys.exit()
            
        someDepth0 = 0
        someLength0 = 0
        for s in segments :
            if s.depth == 0 and s.length > 0:
                if reliableCoverage :
                    if someDepth0 < 10 :
                        print("WARNING: contig ", s.names, " has no readable coverage information or coverage=0. If this is a widespread issue, please use --conservative mode")
                    elif someDepth0 == 10 :
                        print("Not displaying all contigs with no coverage information, but there are more.")
                someDepth0 += 1
            if s.length == 0 :
                s.length1()
                # print("WARNING: contig ", s.names, " has length = 0. This might infer in handling the coverage")
            
        if someDepth0 == len(segments) and reliableCoverage :
            print("WARNING: could not read coverage information in the input GFA. Coverage information for each contig is highly recommended. Continuing nevertheless, switching to --conservative mode")
            reliableCoverage = False
        elif someDepth0>0 and reliableCoverage :
            print("WARNING: ", someDepth0, " contigs out of ", len(segments), " had no coverage information or coverage=0. If this is a widespread issue, please use --conservative mode")

        
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
            
            print("ERROR: You should provide to unzip long reads mapped in GAF format and/or interaction matrices, using either --HiCinteractions (-i) or --linkedReadsInteractions (-k). If you do not have them, you can create them using the HiC-IM or linked-reads-IM commands")
            sys.exit()

        print("================\n\nEverything loaded, moving on to untangling the graph\n\n================")

        #if noisy, clean the graph of small dead-ends and bubbles
        if genomeSize != 0 and reliableCoverage :
            print("Because --genome-size was used, cleaning the graph of small dead-ends and bubbles")
            clean_graph(segments, genomeSize)
            #merge all segments
            merge_adjacent_contigs(segments)
        
        #creating copiesnuber (cn), a dictionnary inventoring how many times each contig appears
        cn = {}
        for segment in segments :
            for name in segment.names :
                cn[name] = 1
        
        ##Moving to the actual unzipping of the graph
        
        supported_links2 = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links considering the topography of the graph
        refHaploidy, multiplicities = determine_multiplicity(segments, supported_links2, reliableCoverage) #multiplicities can be seen as a mininimum multiplicity of each contig regarding the topology of the graph

        #As a first step, use only the long reads, if available
        if uselr :
            print("\n*Untangling the graph using long reads*\n")
            segments = bridge_with_long_reads(segments, names, cn, lrFile, supported_links2, multiplicities, exhaustive)
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)
            print("\n*Done untangling the graph using long reads*\n")
        
        #As a second step, use Hi-C and/or linked reads 
        if interactionMatrix.count_nonzero() > 0 :
            print("\n*Untangling the graph using Hi-C*\n")
            segments = solve_with_HiC(segments, interactionMatrix, names, confidentCoverage=reliableCoverage, noisy = noisy, verbose = verbose, haploid=haploid)
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)
            print("\n*Done untangling the graph using Hi-C*\n")
        
        elif tagInteractionMatrix.count_nonzero() > 0 :
            segments = solve_with_HiC(segments, tagInteractionMatrix, names, confidentCoverage=reliableCoverage, noisy = noisy, verbose = verbose, haploid=haploid)
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)

        elif not uselr :
            
            print("WARNING: all interaction matrices are empty, GraphUnzip does not do anything")
        

            
        # now exporting the output  
        print("Now exporting the result")
        newnames = io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge, rename_contigs=rename)
    
        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile, rename_contigs=rename)
            
        if bamFile != "None":
            print("Now creating the new bam file to re-scaffold")
            io.export_to_bam(segments, bamFile, newnames)
    
        print("Finished in ", time.time() - t, " seconds")
     
    elif command == "purge":
        args = parse_args_purge()
        
        gfaFile = args.gfa
        
        outFile = args.output
        fastaFile = args.fasta_output
        
        merge = not args.dont_merge
        rename = False
        
        segments, names = io.load_gfa(gfaFile)
        
        purge_assembly(segments)
        
        print("Now exporting the result")
        io.export_to_GFA(segments, gfaFile, exportFile=outFile, merge_adjacent_contigs=merge, rename_contigs=rename)
    
        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile, rename_contigs=rename)
            
    
    else :
        print("Unrecognized command ", command, "\". Use either unzip, HiC-IM (to prepare Hi-C data) or linked-reads-IM (to prepare linked reads data)")


if __name__ == "__main__":
    main()
