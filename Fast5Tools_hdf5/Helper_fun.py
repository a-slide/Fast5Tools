# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
from glob import iglob
from collections import OrderedDict

# Third party imports
import numpy as np
import h5py

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

def access_dir (fn, **kwargs):
    """Check if the directory is readeable
    """
    return os.path.isdir (fn) and os.access (fn, os.R_OK),

def access_file (fn, **kwargs):
    """Check if the directory is readeable
    """
    return os.path.isfile (fn) and os.access (fn, os.R_OK),

def parse_attrs (g, **kwargs):
    """Parse hdf5 attributes in an ordered dict"""
    d = OrderedDict ()
    for i, j in g.attrs.items():
        if type(j) == np.bytes_:
            j = j.decode("utf8")
        d[i] = j
    return d

def write_attrs (g, d, **kwargs):
    """Write dict values in an hdf5 file"""
    for i, j in d.items():
        if type(j) == str:
            j = str.encode(j)
        g.attrs.create (name=i, data=j)

def visit_h5 (fn, group=None,**kwargs):
    """visualise the content of a hdf5 file"""
    def printobj(name, obj):
        print(name, obj)
    with h5py.File(fn, "r") as f:
        if not group:
            f.visititems(printobj)
        else:
            f[group].visititems(printobj)
