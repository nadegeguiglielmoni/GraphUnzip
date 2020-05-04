# -*- coding: latin-1 -*-
"""
This file is for handling the gfa : loading it into python lists, transforming it into fasta, reading it...
"""
import time

# from analyse_HiC import lookForStrings
# print(r[:100000])


def gfa_to_fasta():

    gfa_read = open("data/Assembly.gfa")
    r = gfa_read.read()

    seqs = []
    current = ""
    count = 0
    bases = ["A", "C", "G", "T"]

    t1 = time.time()

    for i in range(len(r)):
        count += 1
        if current == "":
            if (
                r[i] in bases
                and r[i + 1] in bases
                and r[i + 2] in bases
                and r[i + 3] in bases
            ):
                current = r[i]

        else:
            # if i%4 == 0 :
            if r[i] in bases:
                current += r[i]
            else:
                seqs += [current]
                current = ""
            # else :
            #   current += r[i]

        if count % 1000000 == 0:
            print(count / 1000000)

        if count == 123533473:
            # print(current)
            break

    fasta_file = open("data/Assembly.fasta", "w")

    for i in range(len(seqs)):
        name = ">sequence" + str(i)
        fasta_file.write(name + "\n")
        fasta_file.write(seqs[i] + "\n")

    print(len(seqs))
    print(len(r))
    print(time.time() - t1)

    return seqs

# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
def load_gfa(file):

    gfa_read = open(file,'r')
    
    names = []
    for line in gfa_read :
        if line[0] == 'S' :
            l = line.split('\t')
            names  += [l[1]]

    links = [[] for i in range(len(names) * 2)]  # a list of links at one end of a contig
    
    gfa_read = open(file,'r')
    
    for line in gfa_read :
        if line[0] == 'L':
            l = line.split('\t')
            contig1 = l[1]
            contig2 = l[3]

            contig1index = names.index(contig1)
            contig2index = names.index(contig2)
            
            if l[2] == '+' and l[4] == '+' :
                links[contig1index*2+1] += [contig2index*2]
                links[contig2index*2] += [contig1index*2+1]
               
            elif l[2] == '+' and l[4] == '-' :
                links[contig1index*2+1] += [contig2index*2+1]
                links[contig2index*2+1] += [contig1index*2+1]
                
            elif l[2] == '-' and l[4] == '-' :
                links[contig1index*2] += [contig2index*2+1]
                links[contig2index*2+1] += [contig1index*2]
                
            elif l[2] == '-' and l[4] == '+' :
                links[contig1index*2] += [contig2index*2]
                links[contig2index*2] += [contig1index*2]
            else :
                print('There seem to be a problem in the gfa file')

    return links, names


def print_short():

    gfa_read = open("data/Assembly.gfa")
    r = gfa_read.read()

    bases = ["A", "C", "G", "T"]
    count = 0
    s = ""

    for i in range(len(r)):
        if r[i] in bases:
            if count % 5000000 == 0:
                s += r[i]
        else:
            s += r[i]
        count += 1
    print(s)
    return 0


# print_short()


# a function to test that our gfa_to_python worked ok
def check_links(links):
    for i, noeud in enumerate(links):
        for j, neighbor in enumerate(noeud):
            found = False
            for k in links[neighbor]:
                if k == i:
                    found = True
            if not found and neighbor != -1:
                print("Problem in links, a one-end link : ", i, neighbor)
                while 1 :
                    r = 0
                return False
    return True
