#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

from graphunzip.segment import Segment
from graphunzip.transform_gfa import check_segments
import graphunzip.input_output as io
import graphunzip.segment as s

# import matplotlib.pyplot as plt
import numpy as np

from bisect import bisect_left  # to look through sorted lists
from copy import deepcopy
import logging
import os


# this function detects and breaks up long (>length) chimeric contigs
def break_up_chimeras(segments, names, interactionMatrix, length):
    allsegments = []
    allXs = []
    for s, segment in enumerate(segments):
        if segment.length > length:
            interactions = []
            X = []

            for axis in range(1, len(segment.names)):
                interaction = 0

                for nameLeft in segment.names[axis:]:
                    for nameRight in segment.names[:axis]:
                        interaction += interactionMatrix[
                            names[nameLeft], names[nameRight]
                        ]

                interactions += [interaction]

                if axis > 1:
                    X += [X[-1] + segment.lengths[axis - 1] / segment.length]
                else:
                    X = [segment.lengths[0] / segment.length]

            allsegments += [interactions]

            # plt.plot (X, interactions)

            inSlump = False
            localMinimums = []
            for axis in range(1, len(interactions) - 1):
                if interactions[axis] < 0.7 * np.max(
                    interactions[:axis]
                ) and interactions[axis] < 0.7 * np.max(interactions[axis:]):
                    if not inSlump:
                        inSlump = True

                    localMinimums += [axis]

                else:
                    if inSlump:
                        inSlump = False

                        loin = [interactions[i] for i in localMinimums].index(
                            np.min([interactions[i] for i in localMinimums])
                        )

                        logging.info(
                            f"Breaking up contig {segment.names} between {segment.names[localMinimums[loin] - 1]} and {segment.names[localMinimums[loin]]} because it looks like a chimeric contig"
                        )

                        # Now break the contig where it should
                        newSegment1, newSegment2 = segment.break_contig(
                            localMinimums[loin]
                        )
                        segments[s] = newSegment1
                        segments.append(newSegment2)

                        localMinimums = []

    # plt.show()
    return segments


# input : one supercontig to be joined with a neighbor at one end
# output : actualized listOfSegments with the two contigs merged
def merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments):
    if len(segment.links[endOfSegment]) != 1:
        logging.error(
            "ERROR : trying to merge simply two contigs that cannot be merged simply."
        )
        return -1, -1

    neighbor = segment.links[endOfSegment][0]
    endOfSegmentNeighbor = segment.otherEndOfLinks[endOfSegment][0]

    if len(neighbor.links[endOfSegmentNeighbor]) != 1:
        logging.error(
            "ERROR : trying to merge simply two contigs that cannot be merged simply."
        )
        return -1, -1

    if neighbor == segment:  # then do not merge a contig with itself
        return -1, -1

    # add the new segment
    s.merge_two_segments(segment, endOfSegment, neighbor, listOfSegments)

    # delete links going towards the two ex-segments
    otherEnd = 1 - endOfSegment
    otherEndNeighbor = 1 - endOfSegmentNeighbor

    for i, n in enumerate(segment.links[otherEnd]):
        n.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)

    for i, n in enumerate(neighbor.links[otherEndNeighbor]):
        # logging.info('Removing ', neighbor.names, ' from ', n.names, ' and adding the new contig',listOfSegments[-1].names, ' at end ', neighbor.otherEndOfLinks[otherEndNeighbor][i])
        n.remove_end_of_link(
            neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor, otherEndNeighbor
        )

    # delete the ex-segments
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)

    return listOfSegments


# a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):
    goOn = True
    while goOn:
        goOn = False
        for segment in listOfSegments:
            alreadyDidThisOne = False  # if the segment is deleted when looking at its first end, you don't want it to look at its other end, since it does not exist anymore
            for endOfSegment in range(2):
                if not alreadyDidThisOne:
                    if (
                        len(segment.links[endOfSegment]) == 1
                        and len(
                            segment.links[endOfSegment][0].links[
                                segment.otherEndOfLinks[endOfSegment][0]
                            ]
                        )
                        == 1
                    ):  # then merge
                        alreadyDidThisOne = True
                        if segment.ID != segment.links[endOfSegment][0].ID:
                            goOn = True
                            listOfSegments = merge_simply_two_adjacent_contig(
                                segment, endOfSegment, listOfSegments
                            )

    return listOfSegments
