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
        elif args[1] == "-v":
            stderr_print ("{} v{}\n".format(package_name, package_version))
        else:
            raise ValueError ("Error: Invalid command '{}'\n".format(args[1]))

    except ValueError as E:
        stderr_print (E)
        stderr_print ("Usage: Fast5Tools [command] [options]\n")
        stderr_print ("Valid command:\n\t-v\n\tfast5_parse\n\tadd_bam_alignment\n\tadd_nanopolish_eventalign\n")
        stderr_print ("For help on given command, type Fast5Tools [command] -h\n")
        sys.exit()

#~~~~~~~~~~~~~~SUBPROGRAMS~~~~~~~~~~~~~~#
def fast5_parse ():
    # Define parser object
    parser = argparse.ArgumentParser(description="""Go get those Fast5!""")
    parser.prog = "Fast5Tools fast5_parse"
    # Define arguments
    parser.add_argument("subprogram", help="Suprogram name = fast5_parse")
    parser.add_argument('--version', '-v', action='version', version=package_version)
    parser.add_argument("-i", "--fast5_dir", required=True, help="Path to the folder containing Fast5 files")
    parser.add_argument("-o", "--db_file", required=True, help="Path to the output database file")
    parser.add_argument("--basecall_group", default='Basecall_1D_000', help="Name of the analysis group in the fast5 file containing the basecalling information")
    parser.add_argument("--raw_read_num", default=0, type=int, help="Index of the read in the Raw")
    parser.add_argument("--threads", default=4, type=int, help="Total number of threads. Minimum = 4")
    parser.add_argument("--verbose", default=False, type=bool, help="If True will be more chatty")
    # Parse Arguments
    a = parser.parse_args()
    # Run command
    f = Fast5Parse(
        fast5_dir = a.fast5_dir,
        db_file=a.db_file,
        basecall_group = a.basecall_group,
        raw_read_num = a.raw_read_num,
        threads = a.threads,
        verbose = a.verbose)

def add_bam_alignment ():
    # Define parser object
    parser = argparse.ArgumentParser(description="""Add Bam alignmnent to a fast5 shelve database""")
    parser.prog = "Fast5Tools add_bam_alignment"
    # Define arguments
    parser.add_argument("subprogram", help="Suprogram name = add_bam_alignment")
    parser.add_argument('--version', '-v', action='version', version=package_version)
    # Parse Arguments
    a = parser.parse_args()
    # Run command
    print ("Not implemented yet")
    parser.print_help ()

def add_nanopolish_eventalign ():
    # Define parser object
    parser = argparse.ArgumentParser(description="""Add Nanopolish Eventalign data to a fast5 shelve database""")
    parser.prog = "Fast5Tools add_nanopolish_eventalign"
    # Define arguments
    parser.add_argument("subprogram", help="Suprogram name = add_nanopolish_eventalign")
    parser.add_argument('--version', '-v', action='version', version=package_version)
    # Parse Arguments
    a = parser.parse_args()
    # Run command
    print ("Not implemented yet")
    parser.print_help ()
