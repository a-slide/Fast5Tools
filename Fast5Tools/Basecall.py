# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np
import pandas as pd

# Local import
from Fast5Tools.Helper_fun import write_attrs

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Basecall (object):
    """
    Represent and summarize basecalling informations from Albacore
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, fastq, kmers, metadata, **kwargs):
        """
        """
        # Self variables
        self.metadata = metadata
        self.kmers = kmers

        # Extract seq and quality from fastq
        fastq_str = fastq.split("\n")
        self.seq = fastq_str[1]
        qual_str = fastq_str[3]
        self.qual = np.empty (shape=len(qual_str), dtype=np.uint32)
        for i, q in enumerate(qual_str):
            self.qual[i] = ord(q)-33

        # Add extra Metadata
        self.metadata["mean_qual"] = self.qual.mean()
        self.metadata["empty_kmers"] = np.isnan(self.kmers["mean"]).sum()

    def __repr__(self):
        """ Readable description of the object """
        m = "[{}]  ".format(self.__class__.__name__)
        if len(self.seq) > 20:
            seq = "{}...{}".format(self.seq[:10], self.seq[-10:])
        else:
            seq = self.seq
        m +="Seq: {} / Length: {} / Empty kmers: {} / Mean quality: {}".format(
            seq, len(self), self.metadata["empty_kmers"], round(self.metadata["mean_qual"], 2))
        return m

    def __len__ (self):
        return len(self.kmers)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#

    @property
    def to_df (self):
        return pd.DataFrame (self.kmers)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _to_hdf5 (self, grp):
        """Write object into an open h5 group"""
        # Save Metadata
        write_attrs (grp, self.metadata)
        # Save Signal
        grp.create_dataset("kmers", data=self.kmers)
        grp.create_dataset("qual", data=self.qual)
        grp.create_dataset("seq", data=str.encode(self.seq))
