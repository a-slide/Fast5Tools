# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from copy import deepcopy
import os
import sys
from collections import OrderedDict

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
        fast5_fn,
        analyses_group='Basecall_1D_000',
        raw_read_num=0,
        basecall = False,
        metadata = False,
        verbose=False,
        **kwargs):
        """
        Parse a Fast5 file basecalled by albacore 2.0+ with h5py and extract the datasets raw, events and fastq.
        The sequence and quality are extracted from the fastq and the event array is collapsed per contiguous kmers
        * fast5_fn: STR
            Path to a fast5 file basecalled by albacore 2.0+
        * analyses_group: STR (default Basecall_1D_000)
            Name of the basecall analyses group in the fast5 file. If None the no basecall values will be fetched
        * raw_read_num: INT (default 0)
            Index of the raw read values in the raw group in the fast5 file. If None the no raw values will be fetched
        """
        # Option self variables
        self.fast5_fn = fast5_fn
        self.basecall = basecall
        self.metadata = metadata

        # Check args
        if not os.access(self.fast5_fn, os.R_OK):
            raise Fast5Error ("Invalid File")

        # Parse the fast5 file
        if verbose: print("Read Fast5 File")
        with h5py.File(fast5_fn, "r") as f:

            # Get raw values
            try:
                self.raw = list(f['/Raw/Reads'].values())[raw_read_num]['Signal'].value
                self.read_id = list(f['/Raw/Reads'].values())[raw_read_num].attrs["read_id"].decode("utf8")
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Raw Value")

            # Get basecall values
            if self.basecall:
                try:
                    fastq = f['/Analyses/{}/BaseCalled_template/Fastq'.format (analyses_group)].value.decode("utf8").split ("\n")
                    self.seq = fastq[1]
                    self.qual = self._qual_str_to_array (fastq[3])
                    events = f['/Analyses/{}/BaseCalled_template/Events'.format (analyses_group)].value
                    self.kmers = self._events_to_kmers (events=events)
                except (KeyError) as E:
                    raise Fast5Error ("No Basecall Value")

            # Get Metadata values
            if self.metadata:
                try:
                    self.context_tags = OrderedDict()
                    for i, j in f['UniqueGlobalKey']["context_tags"].attrs.items():
                        self.context_tags[i] = j.decode("utf8") if type(j) == np.bytes_ else j
                    self.channel_id = OrderedDict()
                    for i, j in f['UniqueGlobalKey']["channel_id"].attrs.items():
                        self.channel_id[i] = j.decode("utf8") if type(j) == np.bytes_ else j
                    self.tracking_id = OrderedDict()
                    for i, j in f['UniqueGlobalKey']["tracking_id"].attrs.items():
                        self.tracking_id[i] = j.decode("utf8") if type(j) == np.bytes_ else j
                except (KeyError) as E:
                    raise Fast5Error ("No Metadata values")

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.fast5_fn)
        m +="\tRead ID: {}\n".format(self.read_id)
        m +="\tCount Raw signals: {}\n".format(self.n_raw)
        if self.basecall:
            m +="\tSequence: {}...\n".format(self.seq[0:25])
            m +="\tMean Read Qual: {}\n".format(round(self.mean_qual, 2))
            m +="\tCount Kmers: {}\n".format(self.n_kmers)
            m +="\tCount Empty Kmers: {}\n".format(self.n_empty_kmers)
        return (m)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def seq_from_kmers (self):
        if self.basecall:
            s=""
            s += self.kmers [0]["seq"][:2]
            for k in self.kmers:
                s += k["seq"][2]
            s += self.kmers [-1]["seq"][3:]
            return s[::-1]

    @property
    def seq_len (self):
        if self.basecall:
            return len(self.seq)

    @property
    def mean_qual (self):
        if self.basecall:
            return self.qual.mean()

    @property
    def n_kmers (self):
        if self.basecall:
            return len(self.kmers)

    @property
    def n_empty_kmers (self):
        if self.basecall:
            return (self.kmers["len"]==0).sum()

    @property
    def n_raw (self):
        return len(self.raw)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot_raw (self,
        start=None,
        end=None,
        kmer_boundaries=False,
        smoothing_win_size=0,
        zscore_norm=False,
        figsize = (30, 5),
        plot_style="ggplot",
        **kwargs):
        """
        Plot raw signal and kmers boundaries
        * start INT
            If defined the raw plot will start at that value
        * end INT
            If defined the raw plot will end at that value
        * kmer_boundaries BOOL
            If True the start and end position of each kmer will be indicated by vertical lines on the graph
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X
        * zscore_norm: BOOL (default False)
            If True the raw will be normalized using the zscore formula
        """

        all_raw = self.raw
        if zscore_norm:
            all_raw = self._zscore_norm (raw=all_raw)
        if smoothing_win_size:
            all_raw = self._raw_signal_smoothing (raw=all_raw, win_size=smoothing_win_size)

        # Define start and end boundaries + extract data
        if not start:
            start = 0
        if not end:
            end=self.n_raw
        raw = all_raw [start:end]
        x_scale = range(start, end)

        # Plot the raw signal

        with pl.style.context(plot_style):
            fig, ax = pl.subplots (figsize=figsize)
            _ = ax.plot (x_scale, raw, color="gray", linewidth=0.5, zorder=1)

            #Plot Kmer boundaries if required and available
            if self.basecall and kmer_boundaries:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.kmers:
                    if row["len"] and row["end"] >= start and row["start"] <= end:
                        y = np.median (all_raw [row["start"]: row["end"]])
                        _ = ax.hlines (y=y, xmin=row["start"], xmax=row["end"], linewidth=1, color="red", zorder=2)

            # Define title
            title = "Read ID {}    Raw:{:,}".format (self.read_id, len(raw))
            if smoothing_win_size: title+= "   Smooth Win Size:{}".format (smoothing_win_size)
            if zscore_norm: title+= "   Z-score normalised"
            _ = ax.set_title (title)
            _ = ax.set_xlim(start, end)

        return fig, ax

    def get_raw (self, start=None, end=None, smoothing_win_size=0, zscore_norm=False):
        """
        * start INT
            If defined the data will start at that value
        * end INT
            If defined the data will end at that value
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X
        * zscore_norm: BOOL (default False)
            If True the raw will be normalized using the zscore formula
        """
        raw = self.raw
        if zscore_norm:
            raw = self._zscore_norm (raw)
        if smoothing_win_size:
            raw = self._raw_signal_smoothing (raw=raw, win_size=smoothing_win_size)
        raw = raw [start:end]

        return raw

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    @classmethod
    def _zscore_norm (cls, raw, **kwargs):
        """ Normalize raw using the zscore formula
        """
        return (raw-raw.mean())/raw.std()

    @classmethod
    def _raw_signal_smoothing (cls, raw, win_size=3, **kwargs):
        """ Smooth the raw signal using a moving median window.
        """
        # Init vals
        n_raw = len(raw)
        #n_smoothed_raw = ((n_raw-win_size)//win_step)+1
        smoothed_raw = np.empty (dtype=np.float, shape=n_raw)

        # Iterate window by window over the raw value array and compute the median for each
        for i, j in enumerate (np.arange (0, n_raw)):
            smoothed_raw [i] = np.median (raw[j:j+win_size])
        return smoothed_raw

    @classmethod
    def _kmer_rephasing (cls, kmers, n_raw, win_step=2, **kwargs):
        """ Rephase kmers coordinates with signal after signal smoothing
        """
        kmers["start"] = kmers["start"]//win_step
        kmers["end"] = kmers["end"]//win_step
        # Corect kmers with index higher that the new raw number, if the last window doesn't reach the end
        kmers["start"][kmers["start"] >= n_raw] = n_raw-1
        kmers["end"][kmers["end"] >= n_raw] = n_raw-1
        # Update the number of empty kmers
        kmers["empty"][kmers["start"] == kmers["end"]] = True
        return kmers

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

    @classmethod
    def _qual_str_to_array (cls, qual_str):
        """ Convert a sanger encoded quality string into a numpy int array. """
        qual = np.zeros(shape=len(qual_str), dtype=np.uint32)
        for i, q in enumerate(qual_str):
            qual[i] = ord(q)-33
        return qual
