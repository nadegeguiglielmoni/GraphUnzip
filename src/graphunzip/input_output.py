#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""

from graphunzip.segment import compute_copiesNumber
from graphunzip.segment import delete_links_present_twice
from graphunzip.segment import Segment
import graphunzip.pybam.pybam as pybam

from scipy import sparse  # to handle interactionMatrix, which should be sparse
import numpy as np

import copy
import gzip
import logging
import os.path  # to check the existence of files
import pickle  # for writing files and reading them
import re  # to find all numbers in a mixed number/letters string (such as 31M1D4M), to split on several characters (<> in longReads_interactionMatrix)
import shutil  # to remove directories
import sys  # to exit when there is an error and to set recursion limit
import time  # to inform the user on what the programm is doing on a regular basis


# Read bam file of Hi-C aligned on assembly
# Input :
#   bam file, assembly graph
# Output :
#   a matrix of interaction between the ends of each contig
def read_bam(file, names, segments):
    interactionMatrix = sparse.dok_matrix((len(segments) * 2, len(segments) * 2))
    number_of_lines = 0

    name_of_last_read = "nothing"
    last_alignment = []
    for alignment in pybam.read(file):
        # logging.info(alignment.sam)
        ls = alignment.sam.split("\t")
        if ls[0] == name_of_last_read:
            if (
                last_alignment[2] in names
                and ls[2] in names
                and last_alignment[2] != ls[2]
            ):  # not an intra-contig contact
                index1 = names[ls[2]] * 2
                if int(ls[3]) > segments[names[ls[2]]].length / 2:
                    index1 += 1

                index2 = names[last_alignment[2]] * 2
                if (
                    int(last_alignment[3])
                    > segments[names[last_alignment[2]]].length / 2
                ):
                    index2 += 1

                interactionMatrix[index1, index2] += 1
                interactionMatrix[index2, index1] += 1

            name_of_last_read = "nothing"
        else:
            name_of_last_read = ls[0]
            last_alignment = ls
        # logging.info(alignment)

        number_of_lines += 1
        if number_of_lines % 10000 == 0:
            logging.info(f"Processed {number_of_lines} records.")

    # logging.info(interactionMatrix)

    interactionMatrix.tocsr()
    return interactionMatrix


# Export the new bam file, ready to be scaffolded
# Input :
#   old bam file, assembly graph
# Output :
#   a new bam file
def export_to_bam(segments, bamFile, newnames):
    newbam = bamFile.rstrip(".bam") + ".new.bam"
    import pysam

    names = {}
    indices = {}
    head = []

    duplicatedcontigs = set()
    for s in segments:
        for i, n in enumerate(s.names):
            if n in names:
                duplicatedcontigs.add(n)
            else:
                names[n] = (
                    s.orientations[i],
                    np.sum([s.lengths[j] for j in range(i)]),
                    newnames[s.full_name()],
                )
                head.append({"LN": s.length, "SN": newnames[s.full_name()]})
                indices[n] = len(head) - 1

    for n in duplicatedcontigs:
        del names[n]

    oldfile = pysam.AlignmentFile(bamFile, "rb")

    header = {"HD": {"VN": "1.5", "SO": "queryname", "GO": "query"}, "SQ": head}

    # header = { 'HD': {'VN': '1.0'},
    #         'SQ': [{'LN': 1575, 'SN': 'chr1'},{'LN': 1584, 'SN': 'chr2'}] }

    # logging.info(header)

    newfile = pysam.AlignmentFile(newbam, "wb", header=header)
    for read in oldfile.fetch(until_eof=True):
        contig = read.reference_name
        start = int(read.reference_start)
        if contig in names:
            a = pysam.AlignedSegment()
            a.query_name = read.query_name
            a.query_sequence = read.query_sequence
            a.flag = read.flag
            a.reference_id = read.reference_id
            a.mapping_quality = read.mapping_quality
            a.cigar = read.cigar
            a.next_reference_id = read.next_reference_id
            a.next_reference_start = read.next_reference_start
            a.template_length = read.template_length
            a.query_qualities = read.query_qualities
            a.tags = read.tags

            a.reference_id = indices[contig]

            if names[contig][1] == "+":
                a.reference_start = start + int(names[contig][0])
            else:
                a.reference_start = start - int(names[contig][0])
            newfile.write(a)

    oldfile.close()
    newfile.close()


# Read fragments list file
# Input :
#   file : fragments_list.txt
# Output :
#   content : list with contig_id, fragment start, fragment end, fragment length
def read_fragment_list(file, header=True):
    with open(file) as f:
        content = f.readlines()

    if header:
        del content[0]

    # parsing
    # 1: contig_id, 2: fragment_start, 3: fragment_end, 4: fragment_length
    content = [x.strip("\n").split("\t") for x in content]
    res = [[x[1], int(x[2]), int(x[3]), int(x[4])] for x in content]

    return res


def read_info_contig(file):
    with open(file) as f:
        content = f.readlines()
    # parsing
    content = [x.strip("\n").split("\t") for x in content[1:]]
    # 1: contig_id, 2: length, 3: n_frags, 4:cumul_length
    content = [[x[0], int(x[1]), int(x[2]), int(x[3])] for x in content]
    return content


# input : output of hicstuff
# output : a matrix of interaction between the ends of each contig
def interactionMatrix(
    hiccontactsfile, fragmentList, names, segments, header=True
):  # the header refers to the hiccontactsfile
    logging.info("Building the interaction matrix.")
    t = time.time()
    # create interaction matrix of contig vs contig
    # 1 -> [1...N] N contigs
    # ...
    # N -> [1...N]
    interactionMatrix = sparse.dok_matrix((len(segments) * 2, len(segments) * 2))

    with open(hiccontactsfile) as f:
        inFile = f.readlines()

    if header:
        del inFile[0]

    n = 0
    unknowncontacts = 0
    for line in inFile:
        if time.time() - t > 2:
            t = time.time()
            logging.info(f"Built {str(int(n / len(inFile) * 100))}%")

        line = line.strip("\n").split("\t")

        # frag1, frag2, contacts
        contact = [int(line[0]), int(line[1]), int(line[2])]

        # search for contig name corresponding to fragment id
        contig1 = fragmentList[contact[0]][0]
        contig2 = fragmentList[contact[1]][0]

        # search for the index of the contigs in names
        if contig1 in names and contig2 in names:
            if (
                fragmentList[contact[0]][1] <= segments[names[contig1]].length / 2
            ):  # if this fragment is more at the left of the contig
                index1 = names[contig1] * 2
            else:  # if the fragment is more at the right of the contig
                index1 = names[contig1] * 2 + 1

            if (
                fragmentList[contact[1]][1] <= segments[names[contig2]].length / 2
            ):  # if this fragment is more at the left of the contig
                index2 = names[contig2] * 2
            else:  # if the fragment is more at the right of the contig
                index2 = names[contig2] * 2 + 1

            if contig1 != contig2:
                # add contacts to interaction matrix
                interactionMatrix[index1, index2] += contact[2]
                interactionMatrix[index2, index1] += contact[2]

                # adds the HiC coverage to the right contigs
                segments[names[contig1]].HiCcoverage += contact[2]
                segments[names[contig2]].HiCcoverage += contact[2]

        else:
            unknowncontacts += 1

        n += 1

    if unknowncontacts != 0:
        logging.warning(
            f"WARNING: There are {unknowncontacts} out of {n} contacts I did not manage to map : you may want to check if the names of the contigs are consistent throughout your files"
        )

    interactionMatrix.tocsr()

    return interactionMatrix


# input : GAF file (outputted by graphaligner) and parameters telling which line are deemed informative
# output : list of useful lines extracted (['>12>34<2' , '>77<33' ,... ] for example)
def read_GAF(
    gafFile, similarity_threshold, whole_mapping_threshold, lines
):  # a function going through the gaf files and inventoring all useful lines
    gaf = open(gafFile, "r")

    for line in gaf:
        ls = line.strip("\n").split("\t")
        path = ls[
            5
        ]  # in GAF format, the 6th column is the path on which the read matched

        if ls[5].count(">") + ls[5].count("<") > 1:
            if (not "id:f" in ls[-2]) or (
                float(ls[-2].split(":")[-1]) > similarity_threshold
            ):
                if (float(ls[3]) - float(ls[2])) / float(
                    ls[1]
                ) > whole_mapping_threshold:
                    lines += [ls[5]]


# input : TSV file (outputted by SPAligner)
# output : list of sequences of contigs in GAF-like format (['>12>34<2' , '>77<33' ,... ] for example)
def read_TSV(tsv_file, names, lines):
    tsv = open(tsv_file, "r")

    for line in tsv:
        ls = line.strip("\n").split("\t")

        alns = ls[6].split(";")

        for aln in alns:
            contigs = aln.split(",")

            if len(contigs) > 1:
                alignment = ""
                for contig in contigs:
                    if "+" in contig[-1]:
                        alignment += ">"
                    else:
                        alignment += "<"
                    alignment += contig[:-1]

                    if contig[:-1] not in names:
                        logging.error(
                            "ERROR: while reading the .tsv, I am coming across a contig that was not in the .gfa, namely ",
                            contig[:-1],
                            ". I recommend you check that you are using the same GFA that you aligned the long reads on.",
                        )
                        sys.exit()
                lines += [alignment]


def linkedReads_interactionMatrix(sam, names):
    import pysam

    interactionMatrix = sparse.dok_matrix((len(names), len(names)))
    contigsInTag = []
    tags = {}
    numbertag = 0
    samfile = pysam.AlignmentFile(
        sam, "r"
    )  # For bam files need to be indexed (faster access)!
    # NOTE: names can also be found in samfile.reference_names from the .BAM header

    possible_tags = ["BX", "BC"]
    l = 0
    # NOTE: for a bam file, the reads are NOT read in the order they come up in the file but ordered by reference seq
    for record in samfile.fetch():
        # no need to check for header, it is parsed in the AlignmentFile object
        ref_name = record.reference_name

        if ref_name in names:  # that means it matched to a contig in the graph
            contig = names[ref_name]

            barcode = None
            for tag in possible_tags:
                if not record.has_tag(tag):
                    continue
                else:
                    barcode = record.get_tag(tag)
                    # logging.info(barcode)
                    break
            else:  # in case barcode is not found (no break statement reached in for loop)
                if l < 10:
                    logging.info(
                        f"Barcode could not be extracted from record {record}, ignoring..."
                    )
                    l += 1
                if l == 9:
                    logging.info(
                        "Other such lines with unextratable barcodes are present, but I will stop displaying them, I think you get the idea"
                    )
                continue  # continue the for record loop

            if barcode in tags:
                contigsInTag[tags[barcode]].append(contig)
            else:
                tags[tag] = numbertag
                contigsInTag += [[contig]]
                numbertag += 1

    for t in contigsInTag:
        for i in range(len(t)):
            for j in range(len(t)):
                interactionMatrix[t[i], t[j]] += 1

    return interactionMatrix


def load_interactionMatrix(file, listOfSegments, names, HiC=False):
    f = open(file, "rb")
    interactionMatrix = pickle.load(f)

    if interactionMatrix.shape != (len(listOfSegments) * 2, len(listOfSegments) * 2):
        logging.error(
            f"ERROR: the interaction matrix provided ({file}) does not seem to match with the GFA file (different number of contigs). The format of interaction matrices have changed since April 2022: you may want to re-run graphunzip HiC-IM. Exiting"
        )
        sys.exit(1)

    if HiC:
        for segment in listOfSegments:
            for contig in segment.names:
                segment.HiCcoverage += np.sum(interactionMatrix[names[contig]])

    return interactionMatrix


# input : contig ID and fasta file
# output : sequence
def get_contig_FASTA(fastaFile, contig, firstline=0):
    with open(fastaFile) as f:
        lookAtNextLine = False
        linenumber = 0
        for line in f:
            if linenumber >= firstline:
                if lookAtNextLine:
                    return line
                target = ">" + str(contig)
                if target in line:
                    lookAtNextLine = True

            linenumber += 1
    return "In get_contig : the contig you are seeking is not in the fasta file"


# input : contig ID, gfa file and contigOffset, the position of the contig in the GFA file
# output : sequence, and if it is present, the sequencing depth of the contig and the rest of the optional tags that could be present in the input gfa
def get_contig_GFA(gfaFile, contig, contigOffset):
    with open(gfaFile) as f:
        f.seek(contigOffset)
        line = f.readline()
        sline = line.strip("\n").split("\t")

        if len(sline) >= 3 and sline[0] == "S" and (contig in sline[1]):
            extra_tags = ""
            depth = ""
            for f in sline[3:]:
                if "dp" in f or "DP" in f or "KC" in f or "RC" in f:
                    depth = f
                else:
                    extra_tags += f + "\t"

            return sline[2], depth, extra_tags

        else:
            logging.error(
                "ERROR : Problem in the offset file, not pointing to the right lines."
            )

    return "In get_contig : the contig you are seeking is not in the gfa file"


# Input :
#   offset file is for speeding up exportation
#   merge_adjacent_contig is to produce a GFA with contigs merged
def export_to_GFA(
    listOfSegments,
    gfaFile="",
    exportFile="results/newAssembly.gfa",
    offsetsFile="",
    merge_adjacent_contigs=False,
    rename_contigs=False,
):
    newnames = {}
    # compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig
    noOffsets = False
    # logging.info('Offsets : ', offsetsFile)
    if offsetsFile == "":
        noOffsets = True
        offsetsFile = gfaFile.strip(".gfa") + "_offsets.pickle"

    if gfaFile != "" and noOffsets:
        # logging.info("coucou")
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile:
            for line in gfafile:
                sline = line.strip("\n").split("\t")
                if sline[0] == "S":
                    # logging.info('In export_to_GFA : exporting ', sline[1])
                    line_offset[
                        sline[1]
                    ] = offset  # adds pair sline[1]:offset to the dict

                offset += len(line)

        # with open(offsetsFile, 'wb') as o:
        #     pickle.dump(line_offset, o)

    # if gfaFile != '' :
    #     with open(offsetsFile, 'rb') as o:
    #         line_offset = pickle.load(o)

    # logging.info(line_offset)

    # logging.info('Line_offsets computed, launching proper writing of the new GFA')
    # Now that the preliminary work is done, start writing the new gfa file

    # now sort the segments by length, to output at the beginning of the files the longests fragments
    listOfSegments.sort(key=lambda x: x.length, reverse=True)

    f = open(exportFile, "w")

    # compute the copiesnumber
    copies = compute_copiesNumber(listOfSegments)

    # write the sequences and the links within the supercontigs
    t = time.time()

    if merge_adjacent_contigs == False:
        for s, segment in enumerate(listOfSegments):
            if time.time() > t + 1:
                t = time.time()
                logging.info(
                    f"{int(s / len(listOfSegments) * 1000) / 10}% of sequences written"
                )

            for c, contig in enumerate(segment.names):
                f.write("S\t" + contig + "-" + str(segment.copiesnumber[c]) + "\t")
                if gfaFile != "":
                    sequence, depth, extra_tags = get_contig_GFA(
                        gfaFile, contig, line_offset[contig]
                    )
                    # logging.info("Here is the depth I got : ", depth)
                    if depth == "":
                        f.write((sequence + "\t" + extra_tags).rstrip("\t") + "\n")
                    else:
                        newdepth = str(float(depth.split(":")[-1]) / copies[contig])
                        f.write(
                            (
                                sequence
                                + "\t"
                                + ":".join(depth.split(":")[:-1])
                                + ":"
                                + newdepth
                                + "\t"
                                + extra_tags
                            ).rstrip("\t")
                            + "\n"
                        )
                else:
                    f.write("*\n")

                if c > 0:
                    f.write(
                        "L\t"
                        + segment.names[c - 1]
                        + "-"
                        + str(segment.copiesnumber[c - 1])
                    )

                    if segment.orientations[c - 1] == 1:
                        f.write("\t+\t")

                    elif segment.orientations[c - 1] == 0:
                        f.write("\t-\t")

                    f.write(contig + "-" + str(segment.copiesnumber[c]))

                    if segment.orientations[c] == 1:
                        f.write("\t+\t")

                    elif segment.orientations[c] == 0:
                        f.write("\t-\t")

                    f.write(segment.insideCIGARs[c - 1] + "\n")

        logging.info("Done exporting sequences, just a little more time...")
        # then write in the gfa file the links between the ends of supercontigs

        for s, segment in enumerate(listOfSegments):
            if time.time() > t + 1:
                t = time.time()
                logging.info(
                    f"{int(s / len(listOfSegments) * 1000) / 10}% of links written."
                )

            for endOfSegment in range(2):
                for l, neighbor in enumerate(segment.links[endOfSegment]):
                    if (
                        segment.ID <= neighbor.ID
                    ):  # that is to ensure each link is written only once
                        endOfNeighbor = segment.otherEndOfLinks[endOfSegment][l]
                        orientation1, orientation2 = "-", "-"

                        if segment.orientations[-endOfSegment] == endOfSegment:
                            orientation1 = "+"

                        if neighbor.orientations[-endOfNeighbor] != endOfNeighbor:
                            orientation2 = "+"

                        f.write(
                            "L\t"
                            + segment.names[-endOfSegment]
                            + "-"
                            + str(segment.copiesnumber[-endOfSegment])
                            + "\t"
                            + orientation1
                            + "\t"
                            + neighbor.names[-endOfNeighbor]
                            + "-"
                            + str(neighbor.copiesnumber[-endOfNeighbor])
                            + "\t"
                            + orientation2
                            + "\t"
                            + segment.CIGARs[endOfSegment][l]
                            + "\n"
                        )

    # in the case the user prefers having merged contigs as an output
    else:  # if merge_adjacent_contigs == True
        # open a file recording which contigs correspond to which supercontigs (with lines such as supercontig_1 contig_A_contig_B_contig_C). Also store that information in a dictionary
        if rename_contigs:
            splitName = exportFile.split("/")[:-1]
            if len(splitName) > 0:
                fcontigs = open("/".join(splitName) + "supercontigs.txt", "w")
            else:
                fcontigs = open("supercontigs.txt", "w")

            supercontigs = {}
            for s, segment in enumerate(listOfSegments):
                supercontigs[segment.full_name()] = "supercontig_" + str(s)

        for s, segment in enumerate(listOfSegments):
            if time.time() > t + 1:
                t = time.time()
                logging.info(
                    f"{int(s / len(listOfSegments) * 1000) / 10}% of sequences written."
                )

            if rename_contigs:
                f.write(
                    "S\t" + "supercontig_" + str(s) + "\t"
                )  # the name of the contigs are supercontig_i
                newnames[segment.full_name()] = "supercontig_" + str(s)
                fcontigs.write(
                    "supercontig_" + str(s) + "\t" + segment.full_name() + "\n"
                )
            else:
                f.write(
                    "S\t" + segment.full_name() + "\t"
                )  # the name of the contigs are supercontig_i
                newnames[segment.full_name()] = segment.full_name()

            fullDepth = 0

            if gfaFile != "":
                sequence = ""
                for c, contig in enumerate(segment.names):
                    s, depth, extra_tags = get_contig_GFA(
                        gfaFile, contig, line_offset[contig]
                    )
                    if s != "*" and segment.orientations[c] == 0:
                        s = s[::-1]
                        complement_dict = {
                            "A": "T",
                            "C": "G",
                            "T": "A",
                            "G": "C",
                            "*": "*",
                        }
                        s = "".join([complement_dict[base] for base in s])
                    if c > 0:
                        CIGARlength = np.sum(
                            [
                                int(i)
                                for i in re.findall(r"\d+", segment.insideCIGARs[c - 1])
                            ]
                        )

                        s = s[CIGARlength:]
                    if depth != "":
                        fullDepth += (
                            float(depth.split(":")[-1]) / copies[contig]
                        ) * len(s)

                    sequence += s

                if fullDepth == 0:
                    f.write((sequence + "\t" + extra_tags).rstrip("\t") + "\n")
                else:
                    newdepth = str(fullDepth / len(sequence))
                    f.write(sequence + "\tDP:f:" + newdepth + "\n")

            else:
                f.write("*\n")

            for endOfSegment in range(2):
                for n, neighbor in enumerate(segment.links[endOfSegment]):
                    if segment.ID < neighbor.ID:  # to write each link just one
                        orientation1, orientation2 = "+", "+"
                        if endOfSegment == 0:
                            orientation1 = "-"
                        if segment.otherEndOfLinks[endOfSegment][n] == 1:
                            orientation2 = "-"

                        if not rename_contigs:
                            f.write(
                                "L\t"
                                + segment.full_name()
                                + "\t"
                                + orientation1
                                + "\t"
                                + neighbor.full_name()
                                + "\t"
                                + orientation2
                                + "\t"
                                + segment.CIGARs[endOfSegment][n]
                                + "\n"
                            )
                        else:
                            f.write(
                                "L\t"
                                + supercontigs[segment.full_name()]
                                + "\t"
                                + orientation1
                                + "\t"
                                + supercontigs[neighbor.full_name()]
                                + "\t"
                                + orientation2
                                + "\t"
                                + segment.CIGARs[endOfSegment][n]
                                + "\n"
                            )

    return newnames


def export_to_fasta(
    listOfSegments,
    gfaFile,
    exportFile="results/newAssembly.fasta",
    rename_contigs=False,
):
    # compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig

    t = 0
    noOffsets = True
    offsetsFile = gfaFile.strip(".gfa") + "_offsets.pickle"

    if gfaFile != "" and noOffsets:
        # logging.info("coucou")
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile:
            for line in gfafile:
                sline = line.strip("\n").split("\t")
                if sline[0] == "S":
                    # logging.info('In export_to_GFA : exporting ', sline[1])
                    line_offset[
                        sline[1]
                    ] = offset  # adds pair sline[1]:offset to the dict

                offset += len(line)

    # logging.info('Line_offsets computed, launching writing of the fasta')
    # Now that the preliminary work is done, start writing the new fasta file

    f = open(exportFile, "w")

    # compute the copiesnumber
    copies = compute_copiesNumber(listOfSegments)

    # now sort the segments by length, to output at the beginning of the files the longests fragments
    listOfSegments.sort(key=lambda x: x.length, reverse=True)

    # Finally, write the sequences
    for s, segment in enumerate(listOfSegments):
        if time.time() > t + 1:
            t = time.time()
            logging.info(
                f"{int(s / len(listOfSegments) * 1000) / 10}% of sequences written."
            )

        if rename_contigs:
            f.write(">supercontig_" + str(s + 1) + "\n")
        else:
            f.write(">" + segment.full_name() + "\n")

        fullDepth = 0

        sequence = ""
        for c, contig in enumerate(segment.names):
            s, depth, extra_contigs = get_contig_GFA(
                gfaFile, contig, line_offset[contig]
            )
            if s != "*" and segment.orientations[c] == 0:
                s = s[::-1]
                complement_dict = {"A": "T", "C": "G", "T": "A", "G": "C"}
                s = "".join([complement_dict[base] for base in s])
            if c > 0:
                CIGARlength = np.sum(
                    [int(i) for i in re.findall(r"\d+", segment.insideCIGARs[c - 1])]
                )

                s = s[CIGARlength:]
            if depth != "":
                fullDepth += (float(depth.split(":")[-1]) / copies[contig]) * len(s)

            sequence += s

        f.write(sequence + "\n")


# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
# Also returns the list of the contig's names
def load_gfa(file):
    logging.info("Loading contigs")
    gfa_read = open(file, "r")

    segments = []

    index = 0
    names = (
        {}
    )  # names is a dictionary that associates the name of each contig in the gfa with an index (which will correspond later to the one in interactionMatrix and copiesnumber)

    for line in gfa_read:
        if line[0] == "S":
            l = line.strip("\n").split("\t")
            cov = 0

            for element in l:
                if "dp" in element[:2] or "DP" in element[:2] or "rd" in element[:2]:
                    try:
                        cov = float(element.split(":")[-1])
                    except:
                        pass

                elif "RC" in element[:2] or "KC" in element[:2]:
                    try:
                        cov = float(element.split(":")[-1]) / len(l[2])
                    except:
                        pass

            s = Segment([l[1]], [1], [len(l[2])], readCoverage=[cov])
            segments.append(s)
            names[
                s.names[0]
            ] = index  # now this contig (identified by its name) is attached to index
            index += 1

    logging.info("Loading links from: {file}")
    gfa_read = open(file, "r")

    cov = 1
    for line in gfa_read:
        if line[0] == "L":
            l = line.strip("\n").split("\t")

            segments[names[l[1]]].add_link_from_GFA(line, names, segments, 0)
            segments[names[l[3]]].add_link_from_GFA(line, names, segments, 1)

    gfa_read.close()

    # delete_links_present_twice(segments)

    return segments, names
