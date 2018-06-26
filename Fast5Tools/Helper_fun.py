# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
from glob import iglob

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#
def find_fast5_files_gen (fast5_dir, **kwargs):
    """
    Generator returning fast5 files found recursively starting from a given folder.
    The recursivity stops as soon as a file matching the extension is found.
    * fast5_dir: STR
        Path to the folder containing Fast5 files (can be in multiple subfolder)
    """
    # In the case where the folder is a file
    if os.path.isdir(fast5_dir):

        # If matching files in the folder
        fast5_found=False
        for fast5 in iglob (os.path.join(fast5_dir, "*.fast5")):
            yield fast5
            fast5_found=True

        # If no matching file go deeper until a leaf containing fast5 is found
        if not fast5_found:
            for item in os.listdir(fast5_dir):
                for fast5 in find_fast5_files_gen (os.path.join(fast5_dir, item)):
                    yield fast5

def stdout_print (*args):
    """reproduce print with stdout.write
    """
    sys.stdout.write(" ".join(str(a) for a in args))
    sys.stdout.flush()
