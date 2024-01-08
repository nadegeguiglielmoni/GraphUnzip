#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""
from graphunzip.arg_helpers import (
    parse_args_command,
    parse_args_purge,
    parse_args_extract,
    parse_args_unzip,
    parse_args_linked,
    parse_args_HiC,
)
import graphunzip.input_output as io

# import analyse_HiC
from graphunzip.finish_untangling import merge_adjacent_contigs
from graphunzip.solve_with_long_reads import bridge_with_long_reads
from graphunzip.transform_gfa import gfa_to_fasta

# from solve_with_long_reads2 import bridge_with_long_reads2
from graphunzip.determine_multiplicity import determine_multiplicity
from graphunzip.solve_with_HiC import solve_with_HiC

# from segment import check_if_all_links_are_sorted
from graphunzip.purge import purge_assembly

from scipy import sparse
import argparse
import logging
import numpy as np
import os.path
import pickle  # reading and writing files
import sys
import time


def main():
    args_command = parse_args_command()
    command = args_command.command

    # if len(sys.argv) < 1 :
    #     sys.exit()

    t = time.time()

    if command == "HiC-IM":
        args = parse_args_HiC()

        matrixFile = args.matrix
        fragmentsFile = args.fragments
        bamfile = args.bam

        if bamfile == "Empty" and (matrixFile == "Empty" or fragmentsFile == "Empty"):
            print(
                "ERROR: you must provide as input either (a bam file) or (an abs_fragments_weighted.txt and a fragment_list.txt files from hicstuff)"
            )
            sys.exit()

        outputIMH = args.HiC_IM

        gfaFile = args.gfa_graph
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0:
            print("ERROR: could not read the GFA")
            sys.exit()

        if bamfile == "Empty":
            if os.path.exists(fragmentsFile) and os.path.exists(matrixFile):
                fragmentList = io.read_fragment_list(fragmentsFile)

                # Now computing the interaction matrix

                interactionMatrix = io.interactionMatrix(
                    matrixFile, fragmentList, names, segments
                )
                useHiC = True

                # exporting it as to never have to do it again

                print("Exporting Hi-C interaction matrix as ", outputIMH)
                with open(outputIMH, "wb") as o:
                    pickle.dump(interactionMatrix, o)

            else:
                print(
                    "Error: could not find fragments file {0}.".format(fragmentsFile),
                    " or matrix file {0}".format(matrixFile),
                )
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

    elif command == "linked-reads-IM":
        args = parse_args_linked()

        barcodedSAM = args.barcoded_SAM
        outputIMT = args.linked_reads_IM

        gfaFile = args.gfa_graph
        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0:
            print("ERROR: could not read the GFA")
            sys.exit()

        if not os.path.exists(barcodedSAM):
            print("Error: could not find the SAM file.")
            sys.exit(1)

        tagInteractionMatrix = io.linkedReads_interactionMatrix(barcodedSAM, names)

        print("Exporting barcoded interaction matrix as ", outputIMT)
        with open(outputIMT, "wb") as o:
            pickle.dump(tagInteractionMatrix, o)

    elif command == "extract":
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
        if len(segments) == 0:
            print("ERROR: could not read the GFA")
            sys.exit()

        # creating copiesnuber (cn), a dictionnary inventoring how many times each contig appears
        cn = {}
        for segment in segments:
            for name in segment.names:
                cn[name] = 1

        print(
            "================\n\nEverything loaded, moving on to untangling the graph\n\n================"
        )

        supported_links2 = sparse.lil_matrix(
            (len(names) * 2, len(names) * 2)
        )  # supported links considering the topography of the graph
        refHaploidy, multiplicities = determine_multiplicity(
            segments, supported_links2, reliable_coverage=False
        )  # multiplicities can be seen as a mininimum multiplicity of each contig regarding the topology of the graph

        segments = bridge_with_long_reads(
            segments,
            names,
            cn,
            lrFile,
            supported_links2,
            multiplicities,
            exhaustive=True,
            extract=True,
        )
        print("Merging contigs that can be merged...")
        merge_adjacent_contigs(segments)
        print("\n*Done extracting the genome*\n")

        # now exporting the output
        print("Now exporting the result")
        io.export_to_GFA(
            segments,
            gfaFile,
            exportFile=outFile,
            merge_adjacent_contigs=merge,
            rename_contigs=rename,
        )

        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile, rename_contigs=rename)

        print("Finished in ", time.time() - t, " seconds")

    elif command == "unzip":
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

        # clean = args.clean

        # Loading the data
        print("Loading the GFA file")
        segments, names = io.load_gfa(
            gfaFile
        )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
        if len(segments) == 0:
            print("ERROR: could not read the GFA")
            sys.exit()

        someDepth0 = 0
        someLength0 = 0
        for s in segments:
            if s.depth == 0:
                if reliableCoverage:
                    if someDepth0 < 10:
                        print(
                            "WARNING: contig ",
                            s.names,
                            " has no readable coverage information or coverage=0. If this is a widespread issue, please use --conservative mode",
                        )
                    elif someDepth0 == 10:
                        print(
                            "Not displaying all contigs with no coverage information, but there are more."
                        )
                someDepth0 += 1
            if s.length == 0:
                s.length1()
                print(
                    "WARNING: contig ",
                    s.names,
                    " has length = 0. This might infer in handling the coverage",
                )

        if someDepth0 == len(segments) and reliableCoverage:
            print(
                "WARNING: could not read coverage information in the input GFA. Coverage information for each contig is highly recommended. Continuing nevertheless, switching to --conservative mode"
            )
            reliableCoverage = False
        elif someDepth0 > 0 and reliableCoverage:
            print(
                "WARNING: ",
                someDepth0,
                " contigs out of ",
                len(segments),
                " had no coverage information or coverage=0. If this is a widespread issue, please use --conservative mode",
            )

        interactionMatrix = sparse.csr_matrix((len(segments), len(segments)))
        tagInteractionMatrix = sparse.csr_matrix((len(segments), len(segments)))
        useHiC = False
        uselr = False
        useTag = False

        if interactionFileH != "Empty":
            if not os.path.exists(interactionFileH):
                print("ERROR: could not access ", interactionFileH)
                sys.exit(1)

            print("Loading the Hi-C interaction matrix")
            interactionMatrix = io.load_interactionMatrix(
                interactionFileH, segments, names, HiC=True
            )
            useHiC = True

        if lrFile != "Empty":
            if not os.path.exists(lrFile):
                print("ERROR: could not access ", lrFile)
                sys.exit(1)
            uselr = True

        if interactionFileT != "Empty":
            if not os.path.exists(interactionFileT):
                print("ERROR: could not access ", interactionFileT)
                sys.exit(1)

            print("Loading the linked-reads interaction matrix")
            tagInteractionMatrix = io.load_interactionMatrix(
                interactionFileT, segments, names, HiC=False
            )
            useTag = True

        if not (useHiC or uselr or useTag):
            print(
                "ERROR: You should provide to unzip long reads mapped in GAF format and/or interaction matrices, using either --HiCinteractions (-i) or --linkedReadsInteractions (-k). If you do not have them, you can create them using the HiC-IM or linked-reads-IM commands"
            )
            sys.exit()

        print(
            "================\n\nEverything loaded, moving on to untangling the graph\n\n================"
        )

        # creating copiesnuber (cn), a dictionnary inventoring how many times each contig appears
        cn = {}
        for segment in segments:
            for name in segment.names:
                cn[name] = 1

        ##Moving to the actual unzipping of the graph

        supported_links2 = sparse.lil_matrix(
            (len(names) * 2, len(names) * 2)
        )  # supported links considering the topography of the graph
        refHaploidy, multiplicities = determine_multiplicity(
            segments, supported_links2, reliableCoverage
        )  # multiplicities can be seen as a mininimum multiplicity of each contig regarding the topology of the graph

        # As a first step, use only the long reads, if available
        if uselr:
            print("\n*Untangling the graph using long reads*\n")
            segments = bridge_with_long_reads(
                segments,
                names,
                cn,
                lrFile,
                supported_links2,
                multiplicities,
                exhaustive,
            )
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)
            print("\n*Done untangling the graph using long reads*\n")

        # As a second step, use Hi-C and/or linked reads
        if interactionMatrix.count_nonzero() > 0:
            print("\n*Untangling the graph using Hi-C*\n")
            segments = solve_with_HiC(
                segments,
                interactionMatrix,
                names,
                confidentCoverage=reliableCoverage,
                noisy=noisy,
                verbose=verbose,
                haploid=haploid,
            )
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)
            print("\n*Done untangling the graph using Hi-C*\n")

        elif tagInteractionMatrix.count_nonzero() > 0:
            segments = solve_with_HiC(
                segments,
                tagInteractionMatrix,
                names,
                confidentCoverage=reliableCoverage,
                noisy=noisy,
                verbose=verbose,
                haploid=haploid,
            )
            print("Merging contigs that can be merged...")
            merge_adjacent_contigs(segments)

        elif not uselr:
            print(
                "WARNING: all interaction matrices are empty, GraphUnzip does not do anything"
            )

        # now exporting the output
        print("Now exporting the result")
        newnames = io.export_to_GFA(
            segments,
            gfaFile,
            exportFile=outFile,
            merge_adjacent_contigs=merge,
            rename_contigs=rename,
        )

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
        io.export_to_GFA(
            segments,
            gfaFile,
            exportFile=outFile,
            merge_adjacent_contigs=merge,
            rename_contigs=rename,
        )

        if fastaFile != "None":
            io.export_to_fasta(segments, gfaFile, fastaFile, rename_contigs=rename)

    else:
        print(
            "Unrecognized command ",
            command,
            '". Use either unzip, HiC-IM (to prepare Hi-C data) or linked-reads-IM (to prepare linked reads data)',
        )


if __name__ == "__main__":
    main()
