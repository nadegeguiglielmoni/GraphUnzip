import pickle #reading and writing files
import os #to run command lines from python (for blast especially)
from Bio.Blast.Applications import NcbiblastnCommandline #to use blast
from Bio.Blast import NCBIXML
#import subprocess #to run command lines from python


# names containing underscores do not work well. Here is a function to get the '_' out of the gfa
def take_underscores_out_of_gfa_names(fileIn, fileOut) :
    
    fo = open(fileOut, 'w')
    fi = open(fileIn, 'r')
    
    for line in fi :
        ls = line.split('\t')
        if 'S' in ls[0] :
            
            ls[1] = ''.join(ls[1].split('_'))
        fo.write('\t'.join(ls))


def cut_chromosomes(solutionFile, chunks = 2000):
    
    os.system('mkdir cut')
    with open(solutionFile) as sol :
        
        queryfiles = []
        q = open('rien', 'w')
        seq = ''
        bit = 0

        for line in sol :
            if '>' in line :
                seq = ''
                bit = 0
                q.close()
                queryfiles += [line.strip('>').strip('\n')+'_cut.fasta']
                q = open('cut/'+line.strip('>').strip('\n')+'_cut.fasta', 'w')
                print('Cutting ', line.strip('>').strip('\n'))
            else :
                
                char = 0
                l = list(line)
                
                while len(l)-char+len(seq) > chunks :
                    
                    while len(seq) < chunks :
                        seq += l[char]
                        char += 1
                        
                    #print('s ' + seq)
                    q.write('>'+str(bit)+'\n')
                    q.write(seq + '\n')
                    bit += 1
                    seq = ''
                
                seq += "".join(list(line)[char:])
                
        os.remove('rien')
        return queryfiles

def assign_a_chromosome_to_each_contig(assemblyFile, queryfiles, fileOut = 'assign.pickle', chunks = 2000):
    
    #here the reference will be the assemblyFile
    os.system('mkdir tmp')
    os.system('makeblastdb -dbtype nucl -parse_seqids -in ' + assemblyFile + ' -out tmp/dtb')
    
    
    assign = {}
    aFile = open(assemblyFile, 'r')
    
    for line in aFile :
        if '>' in line[0] :
            assign[line.strip('>').strip('\n')] = []
    
    #queryfiles = ['2_cut.fasta', '12_cut.fasta']
    
    #map the pieces of chromosome on the assembly file
    
    for query in queryfiles :
        
        contigsHere = []
        print('Looking at chromosome : ', query)
        blastn_cline = NcbiblastnCommandline(query='cut/'+query, db='tmp/dtb', evalue=0.000000001, outfmt=5, out='tmp/blast_sol'+query.strip('.fasta') + '.xml', task="megablast", qcov_hsp_perc = 0.95, num_alignments=5)
        out, err = blastn_cline()
        print('Blast done')
        
        result_handle = open('tmp/blast_sol'+query.strip('.fasta') + '.xml')
        blast_records = NCBIXML.parse(result_handle)
        
        nowhere = 0
        total = 0
        
        for blast_record in blast_records :
            total += 1
            for alignment in blast_record.alignments:
                nbhits = 0
                for hsp in alignment.hsps:
                    if hsp.identities > 0.9*chunks and nbhits == 0:
                        nbhits += 1
                        if query.rstrip('_cut.fasta') not in assign[alignment.title.rstrip('No definition line')] :
                            assign[alignment.title.rstrip('No definition line')] += [query.rstrip('_cut.fasta')]
                
                if nbhits == 0 :
                    nowhere += 1
        
        print('In assign_a_chromosome_to_each_contig, haven\'t managed to map', nowhere/total, ' % of the chunks')

    fo = open(fileOut, 'wb')
    pickle.dump(assign, fo)
    
    os.system('rm -r tmp')
    
    

def check_phasing(assigned, fastaFile) : # the contigs of the fasta file should be merged with '_' between the names
    
    fi = open(fastaFile, 'r')
    
    phasingErrors = 0
    for line in fi :
        
        if '>' in line[0] :
            listOfContigs = line.strip('>').strip('\n').split('-')
            listOfContigs = [''.join(i.split('_')[1:]) for i in listOfContigs]
            #print(listOfContigs)
            
            chromosomes = set()
            allchromosomes = set()
            first = True
            for contig in listOfContigs :
                
                if contig not in assigned :
                    print('There is a problem in the way the contigs are named, possibly if they contain _ or - in their names')
                else :
                    if assigned[contig] != [] :
                        
                        if first :
                            first = False
                            chromosomes = set(assigned[contig])
                            allchromosomes = set(assigned[contig])
                        elif len(chromosomes) != 0 :
                            chromosomes = chromosomes.intersection( set(assigned[contig]))
                            allchromosomes = allchromosomes.union(set(assigned[contig]))
                            if len(chromosomes) == 0 :
                                print ('\nPhasing error detected at contig ', listOfContigs, ', the contigs belong to ', [assigned[i] for i in listOfContigs])
                                phasingErrors += 1
            
            #print('Contig ', listOfContigs, ' is in chromosome ', chromosomes)
    
    phasingSuccesses = 0
    for i in assigned.keys():
        if len(assigned[i]) > 0 :
            phasingSuccesses += 1
    
    print('\n\nTo summarize ', phasingSuccesses, ' contigs contain no phasing errors and ', phasingErrors, ' contain phasing errors')
        
#first cut each chromosome in chunks
# queryfiles = cut_chromosomes('data_A_vaga_HiFi/Flye/tests/A_vaga_12_chromosomes.fasta', chunks = 4000)
# print(queryfiles) #give the value outputted here to query files in the future
# print('Done cutting')

#then assign to each contig of the original assembly which chromosome(s) it belongs to
queryfiles = ['1_cut.fasta', '2_cut.fasta', '3_cut.fasta', '4_cut.fasta', '5_cut.fasta', '6_cut.fasta', '7_cut.fasta', '8_cut.fasta', '9_cut.fasta', '10_cut.fasta', '11_cut.fasta', '12_cut.fasta']
assign_a_chromosome_to_each_contig('A_vaga_article/Nanopore_Ratatosk/avaga.flye_keep-haplotypes_hifi.ont_ratatosk_all.fasta', queryfiles, 'A_vaga_article/Nanopore_Ratatosk/assign.pickle', chunks = 4000)

#if you have the queryfiles and the assign.pickle already you can skip the two above functions.

#check if the new supercontigs are composed of contigs coming from only one chromosome

# f = open('data_A_vaga_HiFi/Flye/assign.pickle', 'rb')
# assigned = pickle.load(f)
# 
# print(assigned)
# print(assigned['521'],
# assigned['881'],
# assigned['854'],
# assigned['776'])
#check_phasing(assigned, 'data_A_Vaga_PacBio/unzipped_merged.fasta')

