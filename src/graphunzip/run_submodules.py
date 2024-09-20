from graphunzip.arg_helpers import (
    parse_args_HiC,
    parse_args_linked,
    parse_args_unzip,
)
import graphunzip.input_output as io

# import analyse_HiC
from graphunzip.finish_untangling import merge_adjacent_contigs
from graphunzip.solve_with_long_reads import bridge_with_long_reads
from graphunzip.transform_gfa import gfa_to_fasta

# from solve_with_long_reads2 import bridge_with_long_reads2
from graphunzip.determine_multiplicity import determine_multiplicity
from graphunzip.solve_with_HiC import solve_with_HiC

from graphunzip.clean_graph import clean_graph

from scipy import sparse

import logging
import os.path
import pickle
import sys


def run_submodule(command):
    if command == "HiC-IM":
        run_prog_hic()
    elif command == "linked-reads-IM":
        run_prog_linked()
    elif command == "unzip":
        run_prog_unzip()


def run_prog_unzip():
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
    logging.info("Loading the GFA file")
    segments, names = io.load_gfa(
        gfaFile
    )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
    if len(segments) == 0:
        logging.error("ERROR: could not read the GFA")
        sys.exit(1)

    someDepth0 = 0
    someLength0 = 0
    for s in segments:
        if s.depth == 0 and s.length > 1:
            if reliableCoverage :
                if someDepth0 < 10 :
                    logging.warning(
                        f"WARNING: contig {s.names} has no readable coverage information or coverage=0. If this is a widespread issue, please use --conservative mode",
                    )
                elif someDepth0 == 10:
                    logging.info(
                        "Not displaying all contigs with no coverage information, but there are more."
                    )
            someDepth0 += 1
        if s.length == 0:
            s.length1()

        if someDepth0 == len(segments) and reliableCoverage :
            print("WARNING: could not read coverage information in the input GFA. Coverage information for each contig is highly recommended. Continuing nevertheless, switching to --conservative mode")
            reliableCoverage = False

        

    if someDepth0 == len(segments) and reliableCoverage:
        logging.warning(
            "WARNING: could not read coverage information in the input GFA. Coverage information for each contig is highly recommended. Continuing nevertheless, switching to --conservative mode"
        )
        reliableCoverage = False
    elif someDepth0 > 0 and reliableCoverage:
        logging.warning(
            "WARNING: {0} contigs out of {1} had no coverage information or coverage=0. If this is a widespread issue, please use --conservative mode".format(someDepth0, len(segments))
        )

    interactionMatrix = sparse.csr_matrix((len(segments), len(segments)))
    tagInteractionMatrix = sparse.csr_matrix((len(segments), len(segments)))
    useHiC = False
    uselr = False
    useTag = False

    if interactionFileH != "Empty":
        if not os.path.exists(interactionFileH):
            logging.error("ERROR: could not access ", interactionFileH)
            sys.exit(1)

        logging.info("Loading the Hi-C interaction matrix")
        interactionMatrix = io.load_interactionMatrix(
            interactionFileH, segments, names, HiC=True
        )
        useHiC = True

    if lrFile != "Empty":
        if not os.path.exists(lrFile):
            logging.error("ERROR: could not access ", lrFile)
            sys.exit(1)
        uselr = True

    if interactionFileT != "Empty":
        if not os.path.exists(interactionFileT):
            logging.error("ERROR: could not access ", interactionFileT)
            sys.exit(1)

        logging.info("Loading the linked-reads interaction matrix")
        tagInteractionMatrix = io.load_interactionMatrix(
            interactionFileT, segments, names, HiC=False
        )
        useTag = True

    if not (useHiC or uselr or useTag):
        logging.error(
            "ERROR: You should provide to unzip long reads mapped in GAF format and/or interaction matrices, using either --HiCinteractions (-i) or --linkedReadsInteractions (-k). If you do not have them, you can create them using the HiC-IM or linked-reads-IM commands"
        )
        sys.exit(1)

    logging.info(
        "================\n\nEverything loaded, moving on to untangling the graph\n\n================"
    )

    #if noisy, clean the graph of small dead-ends and bubbles
    if genomeSize != 0 and reliableCoverage :
        logging.info("Because --genome-size was used, cleaning the graph of small dead-ends and bubbles")
        clean_graph(segments, genomeSize)
        #merge all segments
        merge_adjacent_contigs(segments)

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
        logging.info("Untangling the graph using long reads")
        segments = bridge_with_long_reads(
            segments,
            names,
            cn,
            lrFile,
            supported_links2,
            multiplicities,
            exhaustive,
        )
        logging.info("Merging contigs that can be merged.")
        merge_adjacent_contigs(segments)
        logging.info("Done untangling the graph using long reads.")

    # As a second step, use Hi-C and/or linked reads
    if interactionMatrix.count_nonzero() > 0:
        logging.info("Untangling the graph using Hi-C")
        segments = solve_with_HiC(
            segments,
            interactionMatrix,
            names,
            confidentCoverage=reliableCoverage,
            noisy=noisy,
            verbose=verbose,
            haploid=haploid,
        )
        logging.info("Merging contigs that can be merged...")
        merge_adjacent_contigs(segments)
        logging.info("Done untangling the graph using Hi-C")

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
        logging.info("Merging contigs that can be merged...")
        merge_adjacent_contigs(segments)

    elif not uselr:
        logging.warning(
            "WARNING: All interaction matrices are empty, GraphUnzip does not do anything"
        )

    # now exporting the output
    logging.info("Now exporting the result")
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
        logging.info("Now creating the new bam file to re-scaffold")
        io.export_to_bam(segments, bamFile, newnames)


def run_prog_hic():
    args = parse_args_HiC()

    matrixFile = args.matrix
    fragmentsFile = args.fragments
    bamfile = args.bam

    if bamfile == "Empty" and (matrixFile == "Empty" or fragmentsFile == "Empty"):
        logging.error(
            "ERROR: you must provide as input either (a bam file) or (an abs_fragments_weighted.txt and a fragment_list.txt files from hicstuff)"
        )
        sys.exit(1)

    outputIMH = args.HiC_IM

    gfaFile = args.gfa_graph
    # Loading the data
    logging.info(f"Loading GFA file: {args.gfa_graph}")
    segments, names = io.load_gfa(
        gfaFile
    )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
    if len(segments) == 0:
        logging.error("ERROR: could not read the GFA")
        sys.exit(1)

    if bamfile == "Empty":
        if os.path.exists(fragmentsFile) and os.path.exists(matrixFile):
            fragmentList = io.read_fragment_list(fragmentsFile)

            # Now computing the interaction matrix

            interactionMatrix = io.interactionMatrix(
                matrixFile, fragmentList, names, segments
            )
            useHiC = True

            # exporting it as to never have to do it again

            logging.info("Exporting Hi-C interaction matrix as %s", outputIMH)
            with open(outputIMH, "wb") as o:
                pickle.dump(interactionMatrix, o)

        else:
            logging.error(
                "Error: could not find fragments file {0}.".format(fragmentsFile),
                " or matrix file {0}".format(matrixFile),
            )
            sys.exit(1)
    else:
        if os.path.exists(bamfile):
            interactionMatrix = io.read_bam(bamfile, names, segments)
            useHiC = True

            logging.info("Exporting Hi-C interaction matrix as %s", outputIMH)
            with open(outputIMH, "wb") as o:
                pickle.dump(interactionMatrix, o)
        else:
            logging.error("Error: could not find bam file {0}.".format(bamfile))
            sys.exit(1)


def run_prog_linked():
    args = parse_args_linked()

    barcodedSAM = args.barcoded_SAM
    outputIMT = args.linked_reads_IM

    gfaFile = args.gfa_graph
    # Loading the data
    logging.info(f"Loading GFA file: {args.gfa_graph}")
    segments, names = io.load_gfa(
        gfaFile
    )  # outputs the list of segments as well as names, which is a dict linking the names of the contigs to their index in interactionMatrix, listOfContigs...
    if len(segments) == 0:
        logging.error("ERROR: could not read the GFA")
        sys.exit(1)

    if not os.path.exists(barcodedSAM):
        logging.error("Error: could not find the SAM file.")
        sys.exit(1)

    tagInteractionMatrix = io.linkedReads_interactionMatrix(barcodedSAM, names)

    logging.info("Exporting barcoded interaction matrix as ", outputIMT)
    with open(outputIMT, "wb") as o:
        pickle.dump(tagInteractionMatrix, o)
