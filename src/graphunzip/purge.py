#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 09:25:27 2022

@author: rfaure
"""

from graphunzip.determine_multiplicity import determine_multiplicity

import logging


def purge_assembly(segments):
    IDs = {}
    for s, seg in enumerate(segments):
        IDs[seg.ID] = s

    confidentCoverage = True

    length_of_uncovered_segments = 0.0
    total_length = 0.0
    for s in segments:
        total_length += s.length
        if s.depth == 1 or s.depth == 0:
            length_of_uncovered_segments += s.length

    if float(length_of_uncovered_segments) / total_length > 0.1:
        confidentCoverage = False

    refCoverage, multiplicities = determine_multiplicity(
        segments, reliable_coverage=confidentCoverage
    )

    # now compute the target_ploidy of the input sequence
    ploidies = [0 for i in range(max(multiplicities) + 1)]
    for m in range(len(multiplicities)):
        ploidies[multiplicities[m]] += segments[m].length

    target_ploidy = len(ploidies) - 1
    while target_ploidy >= 0 and ploidies[target_ploidy] < 0.1 * total_length:
        target_ploidy -= 1

    if target_ploidy == -1:
        target_ploidy = ploidies.index(max(ploidies))

    # go through the graph to look for bubbles
    toDelete = []
    for s, seg in enumerate(segments):
        for end in range(2):
            if (
                multiplicities[s] == target_ploidy
            ):  # and "edge_1498" in seg.names and end == 1 :
                unavoidable_segments = unavoidable_neighbors(seg, end, 0, 10, set())

                # logging.info("Looking at contig ", seg.names, " towards end ", end, " and the touching contigs are : ")
                # logging.info([i[0].names for i in unavoidable_segments ])

                unavoidable_segments = {
                    i
                    for i in unavoidable_segments
                    if multiplicities[IDs[i[0].ID]] == target_ploidy
                }
                # list all paths going from the segment to the next unavoidable segment

                if len(unavoidable_segments) > 0:
                    allpaths = find_all_paths(seg, end, unavoidable_segments, set())

                    pathlengths = [sum([i.length for i in j]) for j in allpaths]
                    longestPath = pathlengths.index(
                        max(pathlengths)
                    )  # keep only the longest path

                    for a, alternativePath in enumerate(allpaths):
                        if a != longestPath:
                            for alt_contig in alternativePath:
                                if alt_contig not in allpaths[longestPath]:
                                    alt_contig.cut_all_links()
                                    logging.info("Deleting ", alt_contig.names[0])
                                    toDelete += [IDs[alt_contig.ID]]

    toDelete.sort(reverse=True)
    for t in toDelete:
        del segments[t]


# input: an end of segment and a max recursion depth
# output: the set of contigs you will necessarily find if
def unavoidable_neighbors(segment, end, depth, maxdepth, path):
    if depth >= maxdepth:
        return {(segment, end)}

    result = set()
    first = True
    for n, neighbor in enumerate(segment.links[end]):
        otherEnd = segment.otherEndOfLinks[end][n]
        if (neighbor, 1 - otherEnd) not in path:
            newpath = path.copy()
            newpath.add((segment, end))
            unavoidablethere = unavoidable_neighbors(
                neighbor, 1 - otherEnd, depth + 1, maxdepth, newpath
            )

            if unavoidablethere != set():
                if first:
                    result = unavoidablethere
                    first = False
                else:
                    result = result.intersection(unavoidablethere)

    if first == True and len(segment.links[end]) > 0:
        return set()

    # logging.info("Result for contig ", segment.names, " : ", [i[0].names for i in result])

    if depth != 0:
        result.add((segment, end))
    return result


def find_all_paths(segment, end, unavoidable_segments, seen_contigs):
    result = []
    seen_contigs.add((segment, end))
    for n, neighbor in enumerate(segment.links[end]):
        otherEnd = segment.otherEndOfLinks[end][n]
        if (neighbor, 1 - otherEnd) not in seen_contigs:
            if (neighbor, 1 - otherEnd) in unavoidable_segments:
                return [{segment}]

            pathsthere = find_all_paths(
                neighbor, 1 - otherEnd, unavoidable_segments, seen_contigs
            )
            for i in pathsthere:
                i.add(segment)
            result += pathsthere

    return result
