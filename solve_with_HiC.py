#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 18:06:07 2021

@author: rfaure

A file summarizing the algorithmic part of untangling the graph with HiC
"""

#input : segments and the interactionMatrix of Hi-C contacts. Optionnaly, a list of haploid contigs obtained from the long reads algorithm.
def solve_with_HiC(segments, interactionMatrix, names, list_of_haploid_contigs = []):
    
    haploidArray = [False for i in range(len(segments))] #an array containing "True" for the haploid contigs
    determine_haploid_contigs
    
    list_of_knots = determine_list_of_knots(segments, haploidContigs)

#input: a graph, in the form of a list of segments and a list of haploid conigs
def determine_list_of_knots(segments) :
    return 0