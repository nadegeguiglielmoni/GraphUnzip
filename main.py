#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

@author: zaltabar
"""

import basic_functions as bf
import analyse_HiC
from transform_gfa import load_gfa
from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities

#Loading the data
links, names = load_gfa('data/Assembly.gfa')

"""  #Uncomment if it's the first time the data is ran

fragmentList = bf.read_fragment_list('data/results/fragments_list.txt')

#Now computing the interaction matrix
interactionMatrix = analyse_HiC.interactionMatrix('data/results/abs_fragments_contacts_weighted.txt', fragmentList)

#exporting it as to never have to do it again
bf.export_to_csv(interactionMatrix, 'listsPython/interactionMatrix.csv')

"""

""" #Uncomment if the GFA file does not have a fasta equivalent yet, and if you want to export the output in a GFA complete with sequences

gfa_to_fasta('data/Assembly.gfa','data/Assembly.fasta')

"""

#Importing interactionMatrix if it was already computed
interactionMatrix = bf.import_from_csv('listsPython/interactionMatrix.csv')

stringenceReject = 0.2 #when comparing two links, if the intensity of one is below stringenceReject*intensity of the other, it is deleted
stringenceAccept = 0.45 #when comparing two links, if the intensity of one is above stringenceAccept*intensity of the other, it is confirmed
steps = 10 #number of cycles get_rid_of_bad_links - merge_contigs

links, listOfSuperContigs, copiesnumber = solve_ambiguities(links, interactionMatrix, names, stringenceReject, stringenceAccept, steps)

#now exporting the output
bf.export_to_GFA(links, listOfSuperContigs, copiesnumber, names, 'data/Assempbly.fasta', exportFile = 'results/output.gfa')