# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from copy import deepcopy
import os
import sys

# Third party imports
import numpy as np
import h5py
from matplotlib import pyplot as pl

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
        analyses_group='Basecall_1D_000',
        raw_read_num=0,
        min_read_qual=None,
        min_len=None,
        max_len=None,
        kmer_len=5,
        basecall_required=True,
        smooth_raw_signal = False,
        smoothing_win_size = 10,
        smoothing_win_step = 5,
        verbose=False,
        **kwargs):
        """
        Parse a Fast5 file basecalled by albacore 2.0+ with h5py and extract the datasets raw, events and fastq.
        The sequence and quality are extracted from the fastq and the event array is collapsed per contiguous kmers
        * fast5_file: STR
            Path to a fast5 file basecalled by albacore 2.0+
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
        * basecall_required: BOOL (default True)
            if True will raise an error if no basecall value found
        * smooth_raw_signal: BOOL (default False)
            If True the raw signal will be smoothed and shrink using a moving median window
        * smoothing_win_size: INT (default 10)
            Length of the window used to smooth the raw signal
        * smoothing_win_step: INT (default 5)
            Step of the window used to smooth the raw signal
        """
        # Option self variables
        self.fast5_file = fast5_file

        # Private self variable
        self._smoothing_win_size = smoothing_win_size
        self._smoothing_win_step = smoothing_win_step
        self._kmer_len = kmer_len
        self._smooth_raw_signal = smooth_raw_signal

        # Check args
        if not os.access(self.fast5_file, os.R_OK):
            raise Fast5Error ("Invalid File")

        if verbose: print("Read Fast5 File")
        # Parse the fast5 file
        with h5py.File(fast5_file, "r") as f:

            # Get raw values
            try:
                self.raw = list(f['/Raw/Reads'].values())[raw_read_num]['Signal'].value
                if verbose: print("\tFound raw values")
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Raw Value")

            # Get basecall values
            try:
                events = f['/Analyses/{}/BaseCalled_template/Events'.format(analyses_group)].value
                fastq = f['/Analyses/{}/BaseCalled_template/Fastq'.format(analyses_group)].value
                self._basecall_found = True
                if verbose: print("\tFound basecalling values")
            except (KeyError) as E:
                self._basecall_found = False
                if basecall_required:
                    raise Fast5Error ("No Basecall Value")

        # Processing of basecall values
        if self._basecall_found:
            if verbose: print("Process collected basecalling information")

            # Extract info from fastq sequence and check quality and length
            if verbose: print("\tExtract information from fastq sequence")
            fastq_split = fastq.decode("utf8").split("\n")
            self.seq = fastq_split[1]
            self.qual = self.qual_str_to_array (fastq_split[3])
            if min_read_qual and self.qual.mean() < min_read_qual:
                raise Fast5Error ("Low quality")
            if min_len and len(self.seq) < min_len:
                raise Fast5Error ("Short Sequence")
            if max_len and len(self.seq) > max_len:
                raise Fast5Error ("Long Sequence")

            # Collapse events to kmers
            if verbose: print("\tCollapse events per kmers")
            self.kmers = self._events_to_kmers(events=events, kmer_len=self._kmer_len)

        # Smooth raw signal if required
        if smooth_raw_signal:
            if verbose: print("\tSmooth raw signal")
            self._raw_signal_smoothing(win_size=self._smoothing_win_size, win_step=self._smoothing_win_step)

            assert self.kmers["end"][-1]+1 <= self.n_raw, "Error, kmer end cannot be longer than raw list"

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.fast5_file)
        m +="\tCount Raw signals: {}\n".format(self.n_raw)
        if self._basecall_found:
            m +="\tSequence: {}...\n".format(self.seq[0:25])
            m +="\tMean Read Qual: {}\n".format(round(self.mean_qual, 2))
            m +="\tCount Kmers: {}\n".format(self.n_kmers)
            m +="\tCount Empty Kmers: {}\n".format(self.n_empty_kmers)
        return (m)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def seq_from_kmers (self):
        if self._basecall_found:
            s=""
            for k in self.kmers:
                s = k["seq"][0]+s
            return s

    @property
    def seq_len (self):
        if self._basecall_found:
            return len(self.seq)

    @property
    def mean_qual (self):
        if self._basecall_found:
            return self.qual.mean()

    @property
    def n_kmers (self):
        if self._basecall_found:
            return len(self.kmers)

    @property
    def n_empty_kmers (self):
        if self._basecall_found:
            return (self.kmers["empty"]==True).sum()

    @property
    def n_raw (self):
        return len(self.raw)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot_raw (self, kmer_boundaries=True, **kwargs):
        """
        Plot raw signal and kmers boundaries
        """

        # Plot the raw signal
        ax = pl.subplot()
        _ = ax.plot(self.raw, color="gray", linewidth=0.5)
        _ = ax.set_xlim(0, len(self.raw))

        # Plot Kmer boundaries if required and available
        if kmer_boundaries and self._basecall_found:
            ymin, ymax = ax.get_ylim()
            y1, y2 = ymax, ymax+((ymax-ymin)/10)

            for row in self.kmers:
                if not row["empty"]:
                    _ = ax.vlines(x=row["start"], ymin=y1, ymax=y2, linewidth=0.5, color='green')
                    _ = ax.vlines(x=row["end"], ymin=y1, ymax=y2, linewidth=0.5, color='green')

        # Define title
        title = "Raw:{:,}".format (self.n_raw)
        if self._basecall_found:
            title+= "   Kmers:{:,}   Empty kmers:{:,}   Mean Qual:{}".format (
                self.n_kmers, self.n_empty_kmers, round(self.mean_qual, 2))
        if self._smooth_raw_signal:
            title+= "   Smooth Win Size:{}   Smooth Win Step:{:,}".format (
                self._smoothing_win_size, self._smoothing_win_step)
        _ = ax.set_title (title)

        return ax

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _raw_signal_smoothing (self, win_size=5, win_step=2, **kwargs):
        """
        Smooth the raw signal using a moving median window.
        """
        # Init vals
        n_smoothed_raw = ((self.n_raw-win_size)//win_step)+1
        smoothed_raw = np.empty (dtype=np.int32, shape=n_smoothed_raw)

        # Iterate window by window over the raw value array and compute the median for each
        for i, j in enumerate (np.arange (0, self.n_raw-win_size+1, win_step)):
            smoothed_raw [i] = np.median(self.raw[j:j+win_size])
        self.raw = smoothed_raw

        # Rephase the kmers with the signal if the win_step is more than 1
        if self._basecall_found and win_step > 1:
            self.kmers["start"] = self.kmers["start"]//win_step
            self.kmers["end"] = self.kmers["end"]//win_step
            # Corect kmers with index higher that the new raw number, if the last window doesn't reach the end
            # Update the number of empty kmers
            self.kmers["start"][self.kmers["start"] >= self.n_raw] = self.n_raw-1
            self.kmers["end"][self.kmers["end"] >= self.n_raw] = self.n_raw-1
            self.kmers["empty"][self.kmers["start"] == self.kmers["end"]] = True

    def _events_to_kmers (self, events, kmer_len=5, **kwargs):
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
