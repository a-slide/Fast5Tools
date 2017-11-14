# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from copy import deepcopy
import os
import sys
import warnings

# Third party imports
import numpy as np
import h5py

#~~~~~~~~~~~~~~ERROR DEFINITION~~~~~~~~~~~~~~#
class Fast5Error (Exception):
    """ Generic error for the class
    """
    def __init__(self, *args, **kwargs):
        if not args:
            self.val="default error"
        else:
            self.val=args[0]

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Fast5 (object):
    """
    Parse and extract informations from a Fast5 file basecalled by albacore 2.0+
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
        fast5_file,
        raw_required=True,
        basecall_required=True,
        analyses_group='Basecall_1D_000',
        raw_read_num=0,
        min_read_qual=None,
        min_len=None,
        max_len=None,
        kmer_len=5,
        **kwargs):
        """
        Parse a Fast5 file basecalled by albacore 2.0+ with h5py and extract the datasets raw, events and fastq.
        The sequence and quality are extracted from the fastq and the event array is collapsed per contiguous kmers
        * fast5_file: STR
            Path to a fast5 file basecalled by albacore 2.0+
        * raw_required: BOOL (default True)
            if True will raise an error if no raw value found
        * basecall_required: BOOL (default True)
            if True will raise an error if no basecall value found
        * analyses_group: STR (default Basecall_1D_000)
            Name of the basecall analyses group in the fast5 file. If None the no basecall values will be fetched
        * raw_read_num: INT (default 0)
            Index of the raw read values in the raw group in the fast5 file. If None the no raw values will be fetched
        * min_read_qual: INT or FLOAT  (default None)
            Minimal read quality. If lower, raise a Fast5Error Exception
        * min_len: INT (default None)
            Minimal read length. If shorter, raise a Fast5Error Exception
        * max_len: INT (default None)
            Maximal read lenth. If longer, raise a Fast5Error Exception
        * kmer_len: INT (default 5)
            Length of the kmers in the input data
        """
        # Option self variables
        self.fast5_file = fast5_file
        self.kmer_len = kmer_len
        self.raw_required = raw_required
        self.basecall_required = basecall_required

        # Status self variables
        self.raw_found = False
        self.basecall_found = False

        # Check file Readability
        if not os.access(self.fast5_file, os.R_OK):
            raise Fast5Error ("Invalid File")

        # Parse the fast5 file
        with h5py.File(fast5_file, "r") as f:
            # Get raw values
            try:
                self.raw = list(f['/Raw/Reads'].values())[raw_read_num]['Signal'].value
                self.raw_found = True
            except (KeyError, IndexError, TypeError) as E:
                if self.raw_required:
                    raise Fast5Error ("No Raw Value")
            # Get basecall values
            try:
                self.events = f['/Analyses/{}/BaseCalled_template/Events'.format(analyses_group)].value ## No reason to keep in mem, other than for verification
                self.fastq = f['/Analyses/{}/BaseCalled_template/Fastq'.format(analyses_group)].value
                self.basecall_found = True
            except (KeyError) as E:
                if self.basecall_required:
                    raise Fast5Error ("No Basecall Value")

        # Processing of basecall values
        if self.basecall_found:
            # Extract info from fastq sequence
            fastq_split = self.fastq.decode("utf8").split("\n")
            self.seq = fastq_split[1]
            self.qual_str = fastq_split[3] ## No reason to keep in mem, other than for verification
            self.qual = self.qual_str_to_array (self.qual_str)
            # Check quality and length
            if min_read_qual and self.qual.mean() < min_read_qual:
                raise Fast5Error ("Low quality")
            if min_len and len(self.seq) < min_len:
                raise Fast5Error ("Short Sequence")
            if max_len and len(self.seq) > max_len:
                raise Fast5Error ("Long Sequence")
            # Collapse events per kmers and reconstitute missing kmers
            try:
                self.kmers = self.events_to_kmers (events = self.events, kmer_len= self.kmer_len)
            except Exception as E: ## Unsafe = define the possible errors
                raise Fast5Error ("Kmers collapsing Error")

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.fast5_file)
        if self.raw_found:
            m +="\tCount Raw signals: {}\n".format(self.n_raw)
        if self.basecall_found:
            m +="\tSequence: {}...\n".format(self.seq[0:25])
            m +="\tMean Read Qual: {}\n".format(round(self.mean_qual, 2))
            m +="\tCount Events: {}\n".format(self.n_events)
            m +="\tCount Kmers: {}\n".format(self.n_kmers)
            m +="\tCount Empty Kmers: {}\n".format(self.n_empty_kmers)
        return (m)

    #~~~~~~~~~~~~~~PROPERTIES~~~~~~~~~~~~~~#
    @property
    def seq_from_kmers (self):
        if self.basecall_found:
            s=""
            for k in self.kmers:
                s = k["seq"][0]+s
            return s

    @property
    def seq_len (self):
        if self.basecall_found:
            return len(self.seq)

    @property
    def mean_qual (self):
        if self.basecall_found:
            return self.qual.mean()

    @property
    def n_kmers (self):
        if self.basecall_found:
            return len(self.kmers)

    @property
    def n_empty_kmers (self):
        if self.basecall_found:
            return (self.kmers["empty"]==True).sum()

    @property
    def n_raw (self):
        if self.raw_found:
            return len(self.raw)

    @property
    def n_events (self):
        if self.basecall_found:
            return len(self.events)

    #~~~~~~~~~~~~~~CLASS METHODS~~~~~~~~~~~~~~#
    @classmethod
    def events_to_kmers (cls, events, kmer_len=5, **kwargs):
        """
        Iterate over the event dataframe and merge together contiguous events with the same kmer (move 0).
        Missing kmers are infered from the previous and current kmer sequences.
        * events: numpy ndarray
            2D Array containing the events values obtained from a fast5 file
        * kmer_len: INT (default 5)
            Length of kmers in the source data
        """
        # Filter out 0 move events which doesn't add any info
        events = events[events["move"]>0]
        nmoves = len(events)
        # Create an empty nd array to store elements
        nkmer = events['move'].sum() + kmer_len-1
        kmers = np.empty(
            shape=(nkmer,),
            dtype=[('seq','<U{}'.format(kmer_len)), ('start', '<u8'), ('end', '<u8'), ('empty', np.bool_)])
        kmer_index = 0

        # Iterate over the events ndarray
        for i in np.arange (nmoves):
            move = events["move"][i]
            start = events["start"][i]
            seq = events["model_state"][i].decode("utf8")
            end = (events["start"][i]+events["length"][i]) if i == nmoves-1 else events["start"][i+1]

            if move > 1:
                # Generate the missing kmers by combining the previous and the current sequences
                prev_seq =  "#"*kmer_len if i == 0 else events["model_state"][i-1].decode("utf8")
                for j in range (1, move):
                    missing_kmer_seq = prev_seq [j:move] + seq [0:(kmer_len-move+j)]
                    kmers [kmer_index] = (missing_kmer_seq, start, start, True)
                    kmer_index+=1

            # Add new kmer
            kmers[kmer_index] = (seq, start, end, False)
            kmer_index+=1

        # Add final bases without events
        for i in range (1, kmer_len):
            missing_kmer_seq = seq[i:kmer_len]+"#"*i
            kmers[kmer_index] = (missing_kmer_seq, end, end, True)
            kmer_index+=1

        return kmers

    @classmethod
    def qual_str_to_array (cls, qual_str):
        """ Convert a sanger encoded quality string into a numpy int array. """
        qual = np.zeros(shape=len(qual_str), dtype=np.uint32)
        for i, q in enumerate(qual_str):
            qual[i] = ord(q)-33
        return qual
