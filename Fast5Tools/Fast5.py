# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
from collections import OrderedDict

# Third party imports
import numpy as np
import h5py
from matplotlib import pyplot as pl
import matplotlib.patches as mpatches

# Local import
from Fast5Tools.Basecall import Basecall
from Fast5Tools.Helper_fun import stderr_print

#~~~~~~~~~~~~~~ERROR DEFINITION~~~~~~~~~~~~~~#
class Fast5Error (Exception):
    """ Generic error for the class
    """
    def __init__(self, *args, **kwargs):
        if not args:
            self.err_msg="default Fast5Error"
        else:
            self.err_msg=args[0]

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Fast5 (object):
    """
    Parse and extract informations from a Fast5 file basecalled by albacore 2.0+
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self,
        fast5_fn,
        basecall_group='Basecall_1D_000',
        raw_read_num=0,
        error_on_missing_raw=True,
        error_on_missing_metadata=True,
        error_on_missing_basecall=False,
        signal_normalization = 'zscore',
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
        * signal_normalization (default 'zscore')
            Normalization strategy of the raw signal. Can be None or 'zscore'
        """
        # Self variables
        self.fast5_fn = fast5_fn
        self.signal_normalization = signal_normalization
        self.analyses = OrderedDict()
        self.metadata = OrderedDict()

        # Check args
        if not os.access (self.fast5_fn, os.R_OK):
            raise Fast5Error ("Invalid File")

        # Parse the fast5 file
        with h5py.File(fast5_fn, "r") as f:

            # Get raw values
            try:
                self.raw = list(f['/Raw/Reads'].values())[raw_read_num]['Signal'].value
                if self.signal_normalization == "zscore":
                    self.raw_mean = self.raw.mean()
                    self.raw_std = self.raw.std()
                    self.raw = (self.raw-self.raw_mean)/self.raw_std
                self.read_id = list(f['/Raw/Reads'].values())[raw_read_num].attrs["read_id"].decode("utf8")
            except (KeyError, IndexError, TypeError) as E:
                if error_on_missing_raw:
                    raise Fast5Error ("No Raw Value")

            # Get basecall values
            try:
                fastq = f["Analyses"][basecall_group]["BaseCalled_template"]["Fastq"].value.decode("utf8")
                events = f["Analyses"][basecall_group]["BaseCalled_template"]["Events"].value
                basecall_metadata = OrderedDict ()
                for i, j in f["Analyses"][basecall_group].attrs.items():
                    basecall_metadata[i] = j.decode("utf8") if type(j) == np.bytes_ else j
                self.analyses["Albacore_basecalling"] = Basecall (
                    fastq = fastq,
                    kmers =  self._events_to_kmers (events=events),
                    metadata = basecall_metadata)
            except (KeyError) as E:
                if error_on_missing_basecall:
                    raise Fast5Error ("No Basecall Value")

            # Get Metadata values
            try:
                for tag in ["context_tags", "channel_id", "tracking_id"]:
                    self.metadata[tag] = OrderedDict()
                    for i, j in f["UniqueGlobalKey"][tag].attrs.items():
                        self.metadata[tag][i] = j.decode("utf8") if type(j) == np.bytes_ else j
            except (KeyError) as E:
                if error_on_missing_metadata:
                    raise Fast5Error ("No Metadata values")

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.fast5_fn)
        m +="\tRead ID:{}\tRaw signals points:{}".format(self.read_id, self.n_raw)
        if self.signal_normalization == "zscore":
            m +="\tZscore normalised mean:{}".format(self.raw_mean)
        m += "\n"
        for analyses_name, analysis in self.analyses.items():
            m += "\t{}\n".format(analyses_name)
            m += "\t\t{}\n".format(analysis)
        return (m)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def n_raw (self):
        return len(self.raw)

    @property
    def sequencing_kit (self):
        if "context_tags" in self.metadata and "sequencing_kit" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["sequencing_kit"].upper()

    @property
    def flowcell_type (self):
        if "context_tags" in self.metadata and "flowcell_type" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["flowcell_type"].upper()

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot_raw (self, **kwargs):
        stderr_print ("DEPRECATED, use plot() instead")
        sys.exit()

    def plot (self,
        start=None,
        end=None,
        plot_analyses=["Albacore_basecalling", "Nanopolish_eventalign"],
        smoothing_win_size=0,
        figsize = (30, 5),
        plot_style="ggplot",
        raw_linewidth = 0.5,
        raw_color = "gray",
        **kwargs):
        """
        Plot raw signal and kmers boundaries
        * start INT
            If defined the raw plot will start at that value
        * end INT
            If defined the raw plot will end at that value
        * plot_analyses BOOL
            If True the start and end position of each kmer will be indicated by vertical lines on the graph
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X
        """

        # Define start and end boundaries + extract data
        start = start if start else 0
        end = end if end else self.n_raw
        raw = self.get_raw (start=start, end=end, smoothing_win_size=smoothing_win_size)
        x_scale = range(start, end)

        # Plot the raw signal
        with pl.style.context(plot_style):
            fig, ax = pl.subplots (figsize=figsize)
            _ = ax.plot (x_scale, raw, color=raw_color, linewidth=raw_linewidth, zorder=1)

            legend = []
            #Plot Kmer boundaries if required and available
            if "Albacore_basecalling" in plot_analyses and "Albacore_basecalling" in self.analyses:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.analyses["Albacore_basecalling"].kmers:
                    if row["end"] != row["start"] and row["end"] >= start and row["start"] <= end:
                        _ = ax.hlines (y=row["median"], xmin=row["start"], xmax=row["end"], linewidth=2, color="red", zorder=2, alpha=0.75)
                legend.append (mpatches.Patch(color='red', label='Albacore'))

            if "Nanopolish_eventalign" in plot_analyses and "Nanopolish_eventalign" in self.analyses:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.analyses["Nanopolish_eventalign"].kmers:
                    if row["end"] >= start and row["start"] <= end:
                        _ = ax.hlines (y=row["median"], xmin=row["start"], xmax=row["end"], linewidth=2, color="blue", zorder=3, alpha=0.75)
                legend.append (mpatches.Patch(color='blue', label='Nanopolish'))

            # Define title
            title = "Read ID {} / Raw:{:,}".format (self.read_id, len(raw))
            if smoothing_win_size:
                title+= " / Smooth Win Size:{}".format (smoothing_win_size)
            if self.signal_normalization == "zscore":
                 title+= " / Z-score normalised"
            _ = ax.set_title (title)
            _ = ax.set_xlim(start, end)
            if legend:
                _ = ax.legend (handles=legend, frameon=False, ncol=len(legend))
        return fig, ax

    def get_raw (self, start=None, end=None, smoothing_win_size=0):
        """
        * start INT
            If defined the data will start at that value
        * end INT
            If defined the data will end at that value
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X

        """
        raw = self.raw
        if smoothing_win_size:
            raw = self._raw_signal_smoothing (raw=raw, win_size=smoothing_win_size)
        raw = raw [start:end]
        return raw

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _raw_signal_smoothing (self, raw, win_size=3, **kwargs):
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

    def _events_to_kmers (self, events, **kwargs):
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
            ('seq','<U5'),
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
                prev_seq =  "#"*5 if i == 0 else events["model_state"][i-1].decode("utf8")
                for j in range (1, move):
                    missing_kmer_seq = prev_seq [j:move] + seq [0:(5-move+j)]
                    kmers[kmer_index] = (seq, start, end, np.nan, np.nan, np.nan)
                    kmer_index+=1

            # Add new kmer
            rs = self.raw [start:end]
            kmers[kmer_index] = (seq, start, end, np.mean(rs), np.median(rs), np.std(rs))
            kmer_index+=1

        return kmers
