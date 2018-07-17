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

# Local import
from Fast5Tools.Basecall import Basecall

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
        # Self variables
        self.fast5_fn = fast5_fn
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
                    events = events,
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
        m +="\tRead ID: {}\n".format(self.read_id)
        m +="\tCount Raw signals: {}\n".format(self.n_raw)
        for analyses_name, analysis in self.analyses.items():
            m += "\t{}\n".format(analyses_name)
            m += str(analysis)
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
    def plot_raw (self,
        start=None,
        end=None,
        plot_basecall=True,
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
        * plot_basecall BOOL
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
            if "Albacore_basecalling" in self.analyses and plot_basecall:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.analyses["Albacore_basecalling"].kmers:
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
