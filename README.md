# GraphUnzip

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291093.svg)](https://doi.org/10.5281/zenodo.4291093)

Phases an assembly graph using Hi-C data and/or long reads. 

## Why use GraphUnzip ?

`GraphUnzip` phases an uncollapsed assembly graph in GFA format. Its naive approach makes no assumption on the ploidy or the heterozygosity rate of the organism and thus can be used on highly heterozygous genomes.

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
  -f FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        Optional fasta output [default: None]
  --merge MERGE		If you want to merge in the output file the contigs creating a supercontig into one long contig.
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
  -m MATRIX, --matrix MATRIX
                        Sparse Hi-C contact map
  -F FRAGMENTS, --fragments FRAGMENTS
                        Fragments list
  -i INTERACTIONS, --interactions INTERACTIONS
                        File with interactions [default: None]
  -l LONGREADS, --longreads LONGREADS
                        Long reads mapped to the GFA with GraphAligner (GAF
                        format)
  --exhaustive          Removes all links not found in the GAF file
  -M MINIMUM_MATCH, --minimum_match MINIMUM_MATCH
                        Filters out alignments with a minimum match identity <
                        minimum-match [default: 0]
  -w, --whole_match 	Filters out alignments that do not extend over the
                        whole length of the read
  -d, --debug		Debug mode

```

`GraphUnzip` produces an intermediary file, by default interactionMatrix.pickle. If it is the second time you run `GraphUnzip` on the same dataset, specify with -i the path to this file, it will make the program run much faster.

Recommended options are using -w, --exhaustive and -M with a value corresponding to the precision of the reads (roughly 1-error rate): when using highly precise/corrected reads with an expected error rate of 1% you might want to use -M 0.98, while you might want to use -M 0.7 for high-error rate long reads. The default values of -A and -R should be acceptable for a first run, but you might consider tweaking them:

The accepted threshold is the threshold above which a link is considered real (compared with a competing link). If you notice too many contig duplications, increase this threshold.

The rejected threshold is the threshold below which a link is considered non-existent (compared with a competing link). If the outputted assembly graph is too fragmented, lower this threshold.

## Citation

Please cite `GraphUnzip` using the preprint:

[GraphUnzip: unzipping assembly graphs with long reads and Hi-C](https://www.biorxiv.org/content/10.1101/2021.01.29.428779v1) Roland Faure, Nadège Guiglielmoni and Jean-François Flot, bioRxiv (2020).
doi: https://doi.org/10.1101/2021.01.29.428779
