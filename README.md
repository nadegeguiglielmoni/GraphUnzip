# GraphUnzip

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291093.svg)](https://doi.org/10.5281/zenodo.4291093)

Unzips an assembly graph using Hi-C data and/or long reads and/or linked reads. 

## Why use GraphUnzip ?

`GraphUnzip` improves the contiguity of assembly and duplicates collapsed homozygous contigs, aiming at reconstituting an assembly with haplotypes assembled separately. `GraphUnzip` unzips an uncollapsed assembly graph in GFA format. Its naive approach makes no assumption on the ploidy or the heterozygosity rate of the organism and thus can be used on highly heterozygous genomes. 

## Installation

`GraphUnzip` requires python3 with numpy and scipy, you can install them using `pip install`.
To run `GraphUnzip`, clone this repo using `git clone https://github.com/nadegeguiglielmoni/GraphUnzip.git`, and simply run `main.py`

## Usage

### Input

`GraphUnzip` needs two things to work :

1. An assembly graph in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec) 
and
2. Hi-C data : GraphUnzip needs a sparse contact matrix and a fragment list using the [formats outputted by hicstuff](https://github.com/koszullab/hicstuff#File-formats)
and/or 
2. Long reads (mapped to the GFA in the GAF format of [GraphAligner](https://github.com/maickrau/GraphAligner))
and/or
2. Barcoded linked reads mapped to the contigs of the assembly in [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf). Barcodes need to be designated in the SAM by a BX:Z: tag (e.g. BX:Z:AACTTGTCGGTCAT-1) at the end of each line. A possible pipeline to get this file from barcoded reads using BWA would be to: a) convert the assembly from gfa to fasta format ; b) create a bwa index from the fasta file ; c) align the barcoded reads using BWA with option -C to keep the long reads.


To use `GraphUnzip`, you generally need to proceed in two steps :

1. If using Hi-C or linked reads, build interaction matrix(ces) (a matrix quantifying the pairwise interaction between all contigs): for that use the `HiC-IM`, or `linked-reads-IM` command, depending on which type of data you dispose. You will have to specify the files to which these interaction matrices will be written.
2. Use the command `unzip` to unzip the graph using the interaction matrices built beforehand and/or the gaf file if using long reads. This step is usually extremely quick.


### Options
```bash
python main.py --help

usage: main.py [-h] -g GFA [-o OUTPUT] [-f FASTA_OUTPUT] [-A ACCEPTED]
               [-R REJECTED] [-s STEPS] [-m MATRIX] [-F FRAGMENTS]
               [--HiC_IM HIC_IM] [-i HICINTERACTIONS]
               [-k LINKEDREADSINTERACTIONS] [-l LONGREADS] [-e]
               [--linked_reads_IM LINKED_READS_IM]
               [--barcoded_SAM BARCODED_SAM] [-v] [-d DEBUG] [--merge]
               command

positional arguments:
  command               Either unzip, HiC-IM, long-reads-IM or linked-reads-IM

optional arguments:
  -h, --help            show this help message and exit
  -g GFA, --gfa GFA     GFA file to phase

unzip options:
  -o OUTPUT, --output OUTPUT
                        Output GFA [default: output.gfa]
  -f FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        Optional fasta output [default: None]
  -A ACCEPTED, --accepted ACCEPTED
                        Two links that are compared are deemed both true if
                        the weakest of the two, in term of Hi-C contacts, is
                        stronger than this parameter times the strength of the
                        strongest link [default: 0.30]
  -R REJECTED, --rejected REJECTED
                        When two links are compared, the weakest of the two,
                        in term of Hi-C contacts, is considered false and
                        deleted if it is weaker than this parameter times the
                        strength of the strongest links (always smaller than
                        --accepted)[default: 0.15]
  -s STEPS, --steps STEPS
                        Number of cycles get rid of bad links - duplicate
                        contigs. [default: 10]
  -i HICINTERACTIONS, --HiCinteractions HICINTERACTIONS
                        File containing the Hi-C interaction matrix from HiC-
                        IM [default: None]
  -k LINKEDREADSINTERACTIONS, --linkedReadsInteractions LINKEDREADSINTERACTIONS
                        File containing the linked-reads interaction matrix
                        from linked-reads-IM [default: None]
  -l LONGREADS, --longreads LONGREADS
                        Long reads mapped to the GFA with GraphAligner (GAF
                        format), if you have them
  -v, --verbose
  -d DEBUG, --debug DEBUG
                        Activate the debug mode. Parameter: directory to put
                        the logs and the intermediary GFAs.
  --merge               If you want the output to have all possible contigs
                        merged

HiC-IM options:
  -m MATRIX, --matrix MATRIX
                        Sparse Hi-C contact map
  -F FRAGMENTS, --fragments FRAGMENTS
                        Fragments list
  --HiC_IM HIC_IM       Output file for the Hi-C interaction matrix (required)

linked-reads-IM options:
  --linked_reads_IM LINKED_READS_IM
                        Output file for the linked-read interaction matrix
                        (required)
  --barcoded_SAM BARCODED_SAM
                        SAM file of the barcoded reads aligned to the
                        assembly. Barcodes must still be there (use option -C
                        if aligning with BWA) (required)

```

The default values of -A and -R should be acceptable for a first run, but you might consider tweaking them:

The accepted threshold is the threshold above which a link is considered real (compared with a competing link). If you notice too many contig duplications, increase this threshold.

The rejected threshold is the threshold below which a link is considered non-existent (compared with a competing link). If the outputted assembly graph is too fragmented, lower this threshold.

## Citation

Please cite `GraphUnzip` using the preprint:

[GraphUnzip: unzipping assembly graphs with long reads and Hi-C](https://www.biorxiv.org/content/10.1101/2021.01.29.428779v1) Roland Faure, Nadège Guiglielmoni and Jean-François Flot, bioRxiv (2020).
doi: https://doi.org/10.1101/2021.01.29.428779
