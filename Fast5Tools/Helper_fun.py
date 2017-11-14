# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from glob import iglob, glob
from random import shuffle, sample

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

def find_fast5_files_list (fast5_dir, shuffle_files=False, max_files=None, **kwargs):
    """
    Return a list of fast5 files found recursively starting from a given folder.
    The recursivity stops as soon as a file matching the extension is found.
    * fast5_dir: STR
        Path to the folder containing Fast5 files (can be in multiple subfolder)
    * shuffle_files: BOOL (default False)
        Shuffle files before returning the list
    * max_files: INT (default None)
        if set will return n randomly selected files out of the list
    """
    fast5_list = []
    # In the case where the folder is a file
    if os.path.isdir(fast5_dir):
        # If matching files in the folder
        fast5_list.extend (glob (os.path.join(fast5_dir, "*.fast5")))

        # If no matching file go deeper until a leaf containing fast5 is found
        if not fast5_list:
            for item in os.listdir(fast5_dir):
                fast5_list.extend (find_fast5_files_list (os.path.join(fast5_dir, item)))

    # Postprocessing according to options
    if shuffle_files:
        shuffle(fast5_list)
    if max_files:
        fast5_list= sample(fast5_list, max_files)

    return fast5_list
