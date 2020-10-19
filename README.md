# hic2gfa

Phases an assembly graph using Hi-C data. 

## Why use hic2gfa ?

hic2gfa has the specificity of phasing from an uncollapsed assembly graph in GFA format. Its naive approach makes no assumption on the ploidy or the heterozygosity rate of the organism and thus can be used on highly heterozygous genomes.

## Installation

hic2gfa requires numpy and scipy, you can install them using `pip install`.
## Usage

### Input
hic2gfa needs two things to work :

1. An assembly graph in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec) 
and
2. Hi-C data : hic2gfa needs a sparse contact matrix and a fragment list using the [formats outputted by hicstuff](https://github.com/koszullab/hicstuff#File-formats)
or 
3. Long reads (mapped in the GFA in the GAF format of [GraphAligner](https://github.com/maickrau/GraphAligner))

### Options
```bash
python3 main.py --help
usage: main.py [-h] -g GFA [-o OUTPUT] [-fo FASTA_OUTPUT] [-A ACCEPTED]
               [-R REJECTED] [-s STEPS] [-m MATRIX] [-F FRAGMENTS]
               [-i INTERACTIONS] [-lr LONGREADS]

optional arguments:
  -h, --help            show this help message and exit
  -g GFA, --gfa GFA     GFA file to phase
  -o OUTPUT, --output OUTPUT
                        Output GFA [default: output.gfa]
  -fo FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        Optional fasta output [default: None]
  -A ACCEPTED, --accepted ACCEPTED
                        Threshold to accept links. [default: 0.30]
  -R REJECTED, --rejected REJECTED
                        Threshold to reject links. [default: 0.15]
  -s STEPS, --steps STEPS
                        Number of cycles get rid of bad links - duplicate
                        contigs. [default: 10]
  -m MATRIX, --matrix MATRIX
                        Sparse contact map
  -F FRAGMENTS, --fragments FRAGMENTS
                        Fragments list
  -i INTERACTIONS, --interactions INTERACTIONS
                        File with interactions [default: None]
  -lr LONGREADS, --longreads LONGREADS
                        Long reads mapped to the GFA with GraphAligner (GAF format)

```

`hic2gfa` produces an intermediary file, by default interactionMatrix.pickle. If it is the second time you run `hic2gfa` on the same dataset, specify with -i the path to this file, it will make the program run much faster.

The accepted threshold is the threshold above which a link is considered real (compared with a competing link). If you notice too many contig duplications, increase this threshold.

The rejected threshold is the threshold below which a link is considered non-existent (compared with a competing link). If the outputted assembly graph is too fragmented, lower this threshold.
