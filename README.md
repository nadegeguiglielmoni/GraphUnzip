# GraphUnzip

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291093.svg)](https://doi.org/10.5281/zenodo.4291093)

Unzips an assembly graph using Hi-C data and/or long reads. 

## Why use GraphUnzip ?

`GraphUnzip` unzips an uncollapsed assembly graph in GFA format. Its naive approach makes no assumption on the ploidy or the heterozygosity rate of the organism and thus can be used on highly heterozygous genomes.

## Installation

`GraphUnzip` requires numpy and scipy, you can install them using `pip install`.

## Usage

### Input

`GraphUnzip` needs two things to work :

1. An assembly graph in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec) 
and
2. Hi-C data : GraphUnzip needs a sparse contact matrix and a fragment list using the [formats outputted by hicstuff](https://github.com/koszullab/hicstuff#File-formats)
or 
3. Long reads (mapped to the GFA in the GAF format of [GraphAligner](https://github.com/maickrau/GraphAligner))


To use `GraphUnzip`, you need to proceed in two steps :

1. Build interaction matrix(ces) (a matrix quantifying the pairwise interaction between all contigs): for that use the `HiC-IM` or `long-reads-IM` command, depending on which type of data you dispose. You will have to specify the files to which these interaction matrices will be written.
2. Use the command `unzip` to unzip the graph using the interaction matrices built beforehand. This step is usually extremely quick.


### Options
```bash
python3 main.py --help
python main.py -h
usage: main.py [-h] -g GFA [-o OUTPUT] [-f FASTA_OUTPUT] [-A ACCEPTED]
               [-R REJECTED] [-s STEPS] [-m MATRIX] [-F FRAGMENTS]
               [--HiC_IM HIC_IM] [-i HICINTERACTIONS]
               [-j LONGREADSINTERACTIONS] [-l LONGREADS]
               [--long_reads_IM LONG_READS_IM] [-e] [-M MINIMUM_MATCH] [-w]
               [-v] [-d DEBUG] [--merge]
               command

positional arguments:
  command               Either unzip, HiC-IM or long-reads-IM

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
  -j LONGREADSINTERACTIONS, --longReadsInteractions LONGREADSINTERACTIONS
                        File containing the long-reads interaction matrix from
                        long-reads-IM [default: None]
  -e, --exhaustive      
			Removes all links not found in the GAF file
                        (recommended if you have enough reads)
  -v, --verbose
  -d DEBUG, --debug DEBUG
                        Activate the debug mode. Parameter: directory to put
                        the logs and the intermediary GFAs.
  --merge               
			If you want the output to have all possible contigs
                        merged [default: no]

HiC-IM options:
  -m MATRIX, --matrix MATRIX
                        Sparse Hi-C contact map (required)
  -F FRAGMENTS, --fragments FRAGMENTS
                        Fragments list (required)
  --HiC_IM HIC_IM       
			Output file for the Hi-C interaction matrix (required)

long-reads-IM options:
  -l LONGREADS, --longreads LONGREADS
                        Long reads mapped to the GFA with GraphAligner in the GAF
                        format (required)
  --long_reads_IM LONG_READS_IM
                        Output file for the long-read interaction matrix
                        (required)
  -M MINIMUM_MATCH, --minimum_match MINIMUM_MATCH
                        Filters out alignments with a minimum match identity <
                        minimum-match [default: 0]
  -w, --whole_match     
			Filters out alignments that do not extend over the
                        whole length of the read (recommended if you have
                        enough reads)[default: no]


```

Recommended options are using -w, --exhaustive and -M with a value corresponding to the precision of the reads (roughly 1-error rate): when using highly precise/corrected reads with an expected error rate of 1% you might want to use -M 0.98, while you might want to use -M 0.7 for high-error rate long reads. The default values of -A and -R should be acceptable for a first run, but you might consider tweaking them:

The accepted threshold is the threshold above which a link is considered real (compared with a competing link). If you notice too many contig duplications, increase this threshold.

The rejected threshold is the threshold below which a link is considered non-existent (compared with a competing link). If the outputted assembly graph is too fragmented, lower this threshold.

## Citation

Please cite `GraphUnzip` using the preprint:

[GraphUnzip: unzipping assembly graphs with long reads and Hi-C](https://www.biorxiv.org/content/10.1101/2021.01.29.428779v1) Roland Faure, Nadège Guiglielmoni and Jean-François Flot, bioRxiv (2020).
doi: https://doi.org/10.1101/2021.01.29.428779
