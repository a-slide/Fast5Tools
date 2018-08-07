# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np
import pandas as pd

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
        self.fastq = fastq
        self.metadata = metadata
        self.kmers = kmers

    def __repr__(self):
        """ Readable description of the object """
        m = ""
        if self.seq_len > 20:
            seq = "{}...{}".format(self.fastq_seq[:10], self.fastq_seq[-10:])
        else:
            seq = self.fastq_seq
        m +="\t\tSeq: {} / Length: {} / Empty kmers: {} / Mean quality: {}\n".format(
            seq, self.seq_len, self.n_empty_kmers, round(self.mean_qual, 2))
        return m

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def fastq_seq (self):
        return self.fastq.split("\n")[1]

    @property
    def fastq_qual (self):
        qual_str = self.fastq.split("\n")[3]
        qual = np.zeros(shape=len(qual_str), dtype=np.uint32)
        for i, q in enumerate(qual_str):
            qual[i] = ord(q)-33
        return qual

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
    def mean_qual (self):
        return self.fastq_qual.mean()

    @property
    def n_empty_kmers (self):
        return np.isnan(self.kmers["mean"]).sum()

    @property
    def kmers_df (self):
        return pd.DataFrame (self.kmers)
