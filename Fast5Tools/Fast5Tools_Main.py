  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
from Fast5Tools import __version__ as package_version
from Fast5Tools import __name__ as package_name
from Fast5Tools.Fast5Wrapper import Fast5Wrapper
from Fast5Tools.Fast5Parse import Fast5Parse
from Fast5Tools.Helper_fun import stderr_print

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main ():
    # Simple triage function
    try:
        args = sys.argv
        if len(args) == 1:
            raise ValueError ("Error: Missing command\n")
        if args[1] == "fast5_parse":
            fast5_parse ()
        elif args[1] == "add_bam_alignment":
            add_bam_alignment ()
        elif args[1] == "add_nanopolish_eventalign":
            add_nanopolish_eventalign ()
        elif args[1] in ["-v", "--version"]:
            stderr_print ("{} v{}\n".format(package_name, package_version))
        elif args[1] in ["-h", "--help"]:
            raise ValueError ("Fast5Tools help\n")
        else:
            raise ValueError ("Error: Invalid command '{}'\n".format(args[1]))

    except ValueError as E:
        stderr_print (E)
        stderr_print ("Usage: Fast5Tools [command] [options]\n")
        stderr_print ("Valid command:\n\t-v/--version\n\tfast5_parse\n\tadd_bam_alignment\n\tadd_nanopolish_eventalign\n")
        stderr_print ("For help on given command, type Fast5Tools [command] -h\n")
        sys.exit()

#~~~~~~~~~~~~~~SUBPROGRAMS~~~~~~~~~~~~~~#
def fast5_parse ():
    # Define parser object
    parser = argparse.ArgumentParser(description="fast5_parse read all fast5 files in a directory recursively, extract raw signal, metadata and Albacore basecalling values (if available).\
    The fast5 objects generated are stored in a python GDBM database (shelve)")
    parser.prog = "Fast5Tools fast5_parse"
    # Define arguments
    parser.add_argument("subprogram")
    parser.add_argument("-f", "--fast5_dir", required=True, help="Path to the folder containing Fast5 files")
    parser.add_argument("-d", "--db_file", required=True, help="Path to the output database file")
    parser.add_argument("--basecall_group", default='Basecall_1D_000', help="Name of the analysis group in the fast5 file containing the basecalling information (default = 'Basecall_1D_000')")
    parser.add_argument("--raw_read_num", default=0, type=int, help="Index of the read in the Raw (default = 0)")
    parser.add_argument("--signal_normalization", default="zscore", help="Index of the read in the Raw (default = 'zscore')")
    parser.add_argument("--threads", default=4, type=int, help="Total number of threads (default = 4)")
    parser.add_argument("--verbose", default=False, action='store_true', help="If given will be more chatty (default = False)")
    # Parse Arguments
    a = parser.parse_args()

    # Run command
    f = Fast5Parse(
        fast5_dir = a.fast5_dir,
        db_file=a.db_file,
        basecall_group = a.basecall_group,
        raw_read_num = a.raw_read_num,
        signal_normalization=a.signal_normalization,
        threads = a.threads,
        verbose = a.verbose)

def add_bam_alignment ():
    # Define parser object
    parser = argparse.ArgumentParser(description="Add Bam alignmnent to am existing fast5 shelve database.\
    Hits are saved in an Alignment object associated with the Fast5 read object. Secondary, supplementary and unmapped reads are not saved")
    parser.prog = "Fast5Tools add_bam_alignment"
    # Define arguments
    parser.add_argument("subprogram")
    parser.add_argument("-d", "--db_file", required=True, help="Path to the output database file")
    parser.add_argument("-a", "--alignment_fn", required=True, help="Path to a BAM or SAM file containing all aligned reads (does not need to be sorted, filtered or indexed)")
    parser.add_argument("--analysis_name", help="Name of the analysis in the fast5 file database. Overwrite previous analysis if a same name is already in the db (default = 'Alignment')")
    parser.add_argument("--verbose", default=False, action='store_true', help="If given will be more chatty (default = False)")
    # Parse Argumentss
    a = parser.parse_args()

    # Run command
    f = Fast5Wrapper(
        db_file=a.db_file,
        verbose=a.verbose )

    f.add_bam_alignment (
            alignment_fn = a.alignment_fn,
            analysis_name = a.analysis_name,
            verbose=a.verbose)

def add_nanopolish_eventalign ():
    # Define parser object
    parser = argparse.ArgumentParser(description="Add Nanopolish Eventalign data to am existing fast5 shelve database.\
    Resquiggled read information are saved in an Alignment object associated with the Fast5 read object. Does not support secondary alignment resquiggling\
    Nanopolish eventalign has to run from a bam file containing only the primary alignments with the --print-read-names and --signal-index options")
    parser.prog = "Fast5Tools add_nanopolish_eventalign"
    # Define arguments
    parser.add_argument("subprogram")
    parser.add_argument("-d", "--db_file", required=True, help="Path to the output database file")
    parser.add_argument("-e", "--eventalign_fn", required=True, help="Path to a tsv file output generated by nanopolish eventalign")
    parser.add_argument("--analysis_name", help="Name of the analysis in the fast5 file database, Overwrite previous analysis if a same name is already in the db (default = 'Nanopolish_eventalign')")
    parser.add_argument("--verbose", default=False, action='store_true', help="If given will be more chatty (default = False)")
    # Parse Arguments
    a = parser.parse_args()

    # Run command
    f = Fast5Wrapper(
        db_file=a.db_file,
        verbose=a.verbose )

    f.add_nanopolish_eventalign (
            eventalign_fn=a.eventalign_fn,
            analysis_name=a.analysis_name,
            verbose=a.verbose)
