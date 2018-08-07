# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np
import pandas as pd

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Eventalign (object):
    """
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, ref_name, kmers):
        """"""
        self.ref_name = ref_name
        self.kmers = kmers

    def __repr__(self):
        """ Readable description of the object """
        if self.seq_len > 20:
            seq = "{}...{}".format (self.kmer_seq[:10], self.kmer_seq[-10:])
        else:
            seq = self.self.kmer_seq

        m = "\t\tReference: {} / Sequence: {} / Length: {} \n".format (self.ref_name, seq, self.seq_len)
        return m

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#

    @property
    def kmer_seq (self):
        s=""
        s += self.kmers [0]["seq"][:2]
        for k in self.kmers:
            s += k["seq"][2]
        s += self.kmers [-1]["seq"][3:]
        return s[::-1]

    @property
    def seq_len (self):
        return len(self.kmers)

    @property
    def kmers_df (self):
        return pd.DataFrame (self.kmers)
