# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Basecall (object):
    """
    Represent and summarize basecalling informations
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
        name,
        fastq,
        events,
        metadata,
        **kwargs):
        """
        """
        # Self variables
        self.name = name
        self.fastq = fastq
        self.events = events
        self.metadata = metadata
        self.kmers = self._events_to_kmers (events=events)

    def __repr__(self):
        """ Readable description of the object """
        m="\t[{}] {}\n".format(self.__class__.__name__, self.name)
        m +="\t\tSequence: {}\n".format(self.fastq_seq[0:25])
        m +="\t\tQuality: {}\n".format(self.fastq_qual[0:25])
        m +="\t\tSequence length: {}\tEmpty kmers: {}\tMean quality: {}\n".format(
            self.seq_len, self.n_empty_kmers, self.mean_qual)
        return (m)

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
        return (self.kmers["len"]==0).sum()

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    @classmethod
    def _events_to_kmers (cls, events, **kwargs):
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
        kmers = np.empty (shape=(nkmer,), dtype=[('seq','<U5'), ('start', '<u8'), ('end', '<u8'), ('len', '<u8')])
        kmer_index = 0

        # Iterate over the events ndarray
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
                prev_seq =  "#"*5 if i == 0 else events["model_state"][i-1].decode("utf8")
                for j in range (1, move):
                    missing_kmer_seq = prev_seq [j:move] + seq [0:(5-move+j)]
                    kmers [kmer_index] = (missing_kmer_seq, start, start, 0)
                    kmer_index+=1

            # Add new kmer
            kmers[kmer_index] = (seq, start, end, end-start)
            kmer_index+=1

        return kmers
