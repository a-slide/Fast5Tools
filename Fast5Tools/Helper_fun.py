# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
from glob import iglob

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def stderr_print (*args):
    """reproduce print with stderr.write
    """
    sys.stderr.write(" ".join(str(a) for a in args))
    sys.stderr.flush()

def recursive_file_gen (dir, ext, **kwargs):
    """
    """
    # In the case where the folder is a file
    if os.path.isdir(dir):

        # If matching files in the folder
        file_found=False
        for fn in iglob (os.path.join(dir, "*."+ext)):
            yield fn
            file_found=True

        # If no matching file go deeper until a leaf containing fast5 is found
        if not file_found:
            for item in os.listdir(dir):
                for fn in recursive_file_gen (os.path.join(dir, item), ext):
                    yield fn


def access_dir (dir, **kwargs):
    """Check if the directory is readeable
    """
    return os.path.isdir (dir) and os.access (dir, os.R_OK),
