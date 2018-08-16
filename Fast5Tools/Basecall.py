# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np

# Local import
from Fast5Tools.Helper_fun import write_attrs, parse_attrs

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Basecall (object):
    """
    Represent and summarize basecalling informations from Albacore
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, kmers, seq, qual, metadata, **kwargs):
        """
        """
        # Self variables
        self.kmers = kmers
        self.seq = seq
        self.qual = qual
        self.metadata = metadata

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

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _to_db (self, grp):
        """Write object into an open h5 group"""
        # Save Metadata
        write_attrs (grp, self.metadata)
        # Save Signal
        grp.create_dataset("kmers", data=self.kmers)
        grp.create_dataset("qual", data=self.qual)
        grp.create_dataset("seq", data=str.encode(self.seq))

    #~~~~~~~~~~~~~~CLASS METHODS~~~~~~~~~~~~~~#
    @classmethod
    def from_fast5 (cls, grp, raw):
        # Extract metadata
        metadata = parse_attrs (grp)
        # Extract seq and quality from fastq
        seq, qual = cls._parse_fastq (fastq = grp["BaseCalled_template/Fastq"].value.decode("utf8"))
        # Extract kmers information
        kmers = cls._events_to_kmers (events = grp["BaseCalled_template/Events"].value, raw=raw)
        # Add extra Metadata
        metadata["mean_qual"] = qual.mean()
        metadata["empty_kmers"] = np.isnan(kmers["mean"]).sum()

        return Basecall (kmers=kmers, seq=seq, qual=qual, metadata=metadata)

    @classmethod
    def from_db (cls, grp):
        return Basecall (
            kmers = grp.get("kmers").value,
            seq = grp.get("seq").value.decode("utf8"),
            qual = grp.get("qual").value,
            metadata = parse_attrs (grp))

    @classmethod
    def _events_to_kmers (cls, events, raw, **kwargs):
        """
        Iterate over the event dataframe and merge together contiguous events with the same kmer (move 0).
        Missing kmers are infered from the previous and current kmer sequences.
        * events: numpy ndarray
            2D Array containing the events values obtained from a fast5 file
        """
        # Filter out 0 move events which doesn't add any info
        events = events[events["move"]>0]
        nmoves = len(events)

        # Create an empty nd array to store elements
        nkmer = events['move'].sum()
        kmers = np.empty (shape=(nkmer,), dtype=[
            ('seq','S5'),
            ('start', np.uint32),
            ('end', np.uint32),
            ('mean', np.float64),
            ('median', np.float64),
            ('std', np.float64)])

        # Iterate over the events ndarray
        kmer_index = 0

        for i in np.arange (nmoves):
            move = events["move"][i]
            start = events["start"][i]
            seq = events["model_state"][i].decode("utf8")
            if i == nmoves-1:
                end = events["start"][i] + events["length"][i]
            else:
                end = events["start"][i+1]

            if move > 1:
                # Generate the missing kmers by combining the previous and the current sequences
                if i == 0:
                    prev_seq =  "#"*5
                else:
                    prev_seq = events["model_state"][i-1].decode("utf8")

                for j in range (1, move):
                    missing_kmer_seq = prev_seq [j:move] + seq [0:(5-move+j)]
                    kmers[kmer_index] = (missing_kmer_seq, start, end, np.nan, np.nan, np.nan)
                    kmer_index+=1

            # Add new kmer
            sig = raw.get_signal (start, end)
            kmers[kmer_index] = (seq, start, end, np.mean(sig), np.median(sig), np.std(sig))
            kmer_index+=1

        return kmers

    @classmethod
    def _parse_fastq (cls, fastq):
        """Split fastq in seq and quality"""
        fastq_str = fastq.split("\n")
        seq = fastq_str[1]
        qual_str = fastq_str[3]
        qual = np.empty (shape=len(qual_str), dtype=np.uint32)
        for i, q in enumerate(qual_str):
            qual[i] = ord(q)-33

        return seq, qual
