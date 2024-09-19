import argparse
import sys
from graphunzip._version import __version__


def parse_args_command():
    parser = argparse.ArgumentParser(
        description="Unzips an assembly graph using Hi-C data and/or long reads and/or linked reads.",
        prog="graphunzip",
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    parser.add_argument(
        "command",
        choices=["unzip", "HiC-IM", "linked-reads-IM"],
        help="""
        Sub-command must be one of:
        unzip (untangle the GFA file),
        HiC-IM (to prepare Hi-C data) or
        linked-reads-IM (to prepare linked reads data)
        """,
    )

    return parser.parse_args(sys.argv[1:2])


def parse_args_unzip():
    """
    Gets the arguments from the command line.
    """

    parser = argparse.ArgumentParser()
    groupInput = parser.add_argument_group("Input of GraphUnzip")
    groupOutput = parser.add_argument_group("Output of GraphUnzip")
    groupBehavior = parser.add_argument_group("Behavior of GraphUnzip")
    groupOther = parser.add_argument_group("Other options")
    
    
    groupInput.add_argument("-g", "--gfa", required = True, help="""GFA file to phase""")
 
    groupInput.add_argument(
        "-i",
        "--HiCinteractions",
        required=False,
        default="Empty",
        help="""File containing the Hi-C interaction matrix from HiC-IM [optional]""",
    )
    groupInput.add_argument(
        "-k",
        "--linkedReadsInteractions",
        required=False,
        default="Empty",
        help="""File containing the linked-reads interaction matrix from linked-reads-IM [optional]""",
    )
    groupInput.add_argument(
        "-l", "--longreads", required = False, default="Empty", help="""Long reads mapped to the GFA with GraphAligner (GAF format) or SPAligner (TSV format) [optional]"""
    )

    groupInput.add_argument(
        "-s", "--genomeSize", required = False, default="Empty", help="""Full genome size, counting all haplotypes - e.g. 100m or 3g [optional but recommended]"""
    )

    groupOutput.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.gfa",
        help="""Output GFA [default: output.gfa]""",
    )
    groupOutput.add_argument(
        "-f",
        "--fasta_output",
        required=False,
        default="None",
        help="""Optional fasta output [default: None]""",
    )
    groupOutput.add_argument(
        "-b",
        "--bam_file",
        required=False,
        default="None",
        help="""bam file of the Hi-C reads aligned on assembly. GraphUnzip will output bam_file.new.bam corresponding to the new bam file, ready to be used for scaffolding [optional]""",
    )
    
    groupOther.add_argument(
        "-v",
        "--verbose",
        required = False,
        action="store_true",
    )
    groupOther.add_argument(
        "-r",
        "--dont_rename",
        action="store_true",
        help="""Use if you don't want to name the resulting supercontigs with short names but want to keep the names of the original contigs""",
    )
    # groupOther.add_argument(
    #     "-d",
    #     "--debug",
    #     required = False,
    #     default = '',
    #     help="""Activate the debug mode. Parameter: directory to put the logs and the intermediary GFAs.""",
    # )
    groupOther.add_argument(
        "--dont_merge",
        required=False,
        action="store_true",
        help="""If you don't want the output to have all possible contigs merged""",
    )
    
    groupBehavior.add_argument(
        "-H",
        "--haploid",
        action="store_true",
        help="""Use this option if you wish to obtain a collapsed assembly of a multiploid genome.""",
    )
    groupBehavior.add_argument(
        "-c",
        "--conservative",
        action="store_true",
        help="""(Hi-C only) Output very robust contigs. Use this option if the coverage information of the graph is not reliable""",
    )
    
    groupBehavior.add_argument(
        "-B",
        "--bold",
        action="store_true",
        help="""(Hi-C only)[default] Proposes the best untangling it can get (can be misled by approximate coverage information). Use this option if the contig coverage information of the graph can be trusted""",
    )
    groupBehavior.add_argument(
        "-n",
        "--noisy",
        action="store_true",
        help="""(Hi-C only) Use this option if you expect that the assembly may contain artefactual contigs, e.g. when you use the .p_utg.gfa of hifiasm""",
    )
    groupBehavior.add_argument(
        "-e",
        "--exhaustive",
        action="store_true",
        help="""(long reads only) All links not found in the .gaf will be removed""",
    )

    
    return parser.parse_args(sys.argv[2:])


def parse_args_linked():
    """
    Gets the arguments from the command line.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-g",
        "--gfa_graph",
        required=True,
        help="""GFA file that will be untangled (required)""",
    )

    parser.add_argument(
        "-p",
        "--linked_reads_IM",
        required=True,
        help="""Output file for the linked-read interaction matrix (required)""",
    )

    parser.add_argument(
        "-b",
        "--barcoded_SAM",
        required=True,
        help="""SAM file of the barcoded reads aligned to the assembly. Barcodes must still be there (use option -C if aligning with BWA) (required)""",
    )

    return parser.parse_args(sys.argv[2:])


def parse_args_HiC():
    """
    Gets the arguments from the command line.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-g",
        "--gfa_graph",
        required=True,
        help="""GFA file that will be untangled (required)""",
    )
    parser.add_argument(
        "-b",
        "--bam",
        default="Empty",
        required=False,
        help="""Bam file of Hi-C reads aligned on assembly and sorted by name (if using bam format)""",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        default="Empty",
        required=False,
        help="""Sparse Hi-C contact map (if using instaGRAAL format)""",
    )
    parser.add_argument(
        "-F",
        "--fragments",
        default="Empty",
        required=False,
        help="""Fragments list (if using instaGRAAL format)""",
    )
    parser.add_argument(
        "-i",
        "--HiC_IM",
        required=True,
        help="""Output file for the Hi-C interaction matrix (required)""",
    )

    return parser.parse_args(sys.argv[2:])
