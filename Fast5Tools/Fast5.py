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
from Fast5Tools.Raw import Raw
from Fast5Tools.Helper_fun import stderr_print, parse_attrs, write_attrs

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
        basecall_id='Basecall_1D_000',
        basecall_required=False,
        signal_normalization = 'zscore',
        **kwargs):
        """
        Parse a Fast5 file basecalled by albacore 2.0+ with h5py and extract the datasets raw, events and fastq.
        The sequence and quality are extracted from the fastq and the event array is collapsed per contiguous kmers
        * fast5_fn: STR
            Path to a fast5 file basecalled by albacore 2.0+
        * basecall_id: STR (default Basecall_1D_000)
            Name of the basecall analyses group in the fast5 file. If None the no basecall values will be fetched
        * signal_normalization (default 'zscore')
            Normalization strategy of the raw signal. Can be None or 'zscore'
        """
        # Check args
        if not os.access (fast5_fn, os.R_OK):
            raise Fast5Error ("Invalid File")

        # Parse the fast5 file
        with h5py.File(fast5_fn, "r") as f:

            # Get Metadata values
            try:
                self.metadata = OrderedDict ()
                self.metadata["context_tags"] = parse_attrs (f.get("UniqueGlobalKey/context_tags"))
                self.metadata["tracking_id"] = parse_attrs (f.get("UniqueGlobalKey/tracking_id"))
                self.metadata["channel_id"] = parse_attrs (f.get("UniqueGlobalKey/channel_id"))
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Metadata values")

            # Get raw values
            try:
                raw_grp = list(f['/Raw/Reads'].values())[0]
                self.raw = Raw (
                    signal = raw_grp['Signal'].value,
                    metadata = parse_attrs (raw_grp),
                    normalization = signal_normalization)
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Raw Value")

            # Get basecall values if available
            try:
                basecall_group = f["Analyses"][basecall_id]
                self.basecall = Basecall (
                    fastq = basecall_group["BaseCalled_template/Fastq"].value.decode("utf8"),
                    kmers = self._events_to_kmers( basecall_group["BaseCalled_template/Events"].value),
                    metadata = parse_attrs (basecall_group))
                self.has_basecall = True
            except (KeyError, IndexError, TypeError) as E:
                self.has_basecall = False
                if basecall_required:
                    raise Fast5Error ("No Basecall Value")

    def __repr__(self):
        """ Readable description of the object """
        m="[{}]\tRead ID: {}\n".format(self.__class__.__name__, self.read_id)
        # Raw
        m += "\t{}\n".format(self.raw)
        # Basecall if available
        if self.has_basecall:
            m += "\t{}\n".format(self.basecall)
        return (m)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def read_id (self):
        return self.raw.metadata["read_id"]

    @property
    def sequencing_kit (self):
        if "sequencing_kit" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["sequencing_kit"].upper()

    @property
    def flowcell_type (self):
        if "flowcell_type" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["flowcell_type"].upper()

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot_raw (self, **kwargs):
        stderr_print ("DEPRECATED, use plot() instead")
        sys.exit()

    def plot (self,
        start=None,
        end=None,
        plot_analyses=["Basecall"],
        smoothing_win_size=0,
        figsize = (30, 5),
        plot_style="ggplot",
        signal_linewidth = 0.5,
        signal_color = "gray",
        **kwargs):
        """
        Plot signal and kmers boundaries
        * start INT
            If defined the signal plot will start at that value
        * end INT
            If defined the signal plot will end at that value
        * plot_analyses BOOL
            If True the start and end position of each kmer will be indicated by vertical lines on the graph
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X
        """

        # Define start and end boundaries + extract data
        start = start if start else 0
        end = end if end else len(self.raw)
        signal = self.raw.get_signal (start=start, end=end, smoothing_win_size=smoothing_win_size)
        x_scale = range(start, end)

        # Plot the signal
        with pl.style.context(plot_style):
            fig, ax = pl.subplots (figsize=figsize)
            _ = ax.plot (x_scale, signal, color=signal_color, linewidth=signal_linewidth, zorder=1)

            legend = []
            #Plot Kmer boundaries if required and available
            if "Basecall" in plot_analyses and "basecall" in self.__dict__:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.basecall.kmers:
                    if row["end"] != row["start"] and row["end"] >= start and row["start"] <= end:
                        _ = ax.hlines (y=row["median"], xmin=row["start"], xmax=row["end"], linewidth=2, color="red", zorder=2, alpha=0.75)
                legend.append (mpatches.Patch(color='red', label='Basecall'))

            # if "Nanopolish_eventalign" in plot_analyses and "Nanopolish_eventalign" in self.analyses:
            #     ymin, ymax = ax.get_ylim()
            #     y1, y2 = ymax, ymax+((ymax-ymin)/10)
            #
            #     for row in self.analyses["Nanopolish_eventalign"].kmers:
            #         if row["end"] >= start and row["start"] <= end:
            #             _ = ax.hlines (y=row["median"], xmin=row["start"], xmax=row["end"], linewidth=2, color="blue", zorder=3, alpha=0.75)
            #     legend.append (mpatches.Patch(color='blue', label='Nanopolish'))

            # Define title
            title = "Read ID {} / Raw:{:,}".format (self.read_id, len(signal))
            _ = ax.set_title (title)
            _ = ax.set_xlim(start, end)
            if legend:
                _ = ax.legend (handles=legend, frameon=False, ncol=len(legend))
        return fig, ax


    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
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
            sig = self.raw.get_signal (start, end)
            kmers[kmer_index] = (seq, start, end, np.mean(sig), np.median(sig), np.std(sig))
            kmer_index+=1

        return kmers

    def _to_hdf5 (self, grp):
        """Write object into an open h5 group"""

        # Save Metadata
        #md_grp = grp.create_group("metadata")
        write_attrs (grp.create_group("context_tags"), self.metadata["context_tags"])
        write_attrs (grp.create_group("tracking_id"), self.metadata["tracking_id"])
        write_attrs (grp.create_group("channel_id"), self.metadata["channel_id"])

        # Save raw
        self.raw._to_hdf5 (grp.create_group("raw"))

        #Save basecall
        if self.has_basecall:
            self.basecall._to_hdf5 (grp.create_group("basecall"))
