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
from Fast5Tools_hdf5.Basecall import Basecall
from Fast5Tools_hdf5.Raw import Raw
from Fast5Tools_hdf5.Helper_fun import stderr_print, parse_attrs, write_attrs

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
        read_id,
        metadata,
        data,
        **kwargs):
        """
        """
        self.read_id = read_id
        self.metadata = metadata
        self.data = data

    def __repr__(self):
        """ Readable description of the object """
        m="[{}]\tRead ID: {}\n".format(self.__class__.__name__, self.read_id)
        for i, j in self.data.items():
            m += "\t[{}] {}\n".format(i, j)
        return (m)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#

    @property
    def sequencing_kit (self):
        if "sequencing_kit" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["sequencing_kit"].upper()

    @property
    def flowcell_type (self):
        if "flowcell_type" in self.metadata["context_tags"]:
            return self.metadata["context_tags"]["flowcell_type"].upper()

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def plot (self,
        start=None,
        end=None,
        plot_analyses=["basecall"],
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
        end = end if end else len(self.data["raw"])
        signal = self.data["raw"].get_signal (start=start, end=end, smoothing_win_size=smoothing_win_size)
        x_scale = range(start, end)

        # Plot the signal
        with pl.style.context(plot_style):
            fig, ax = pl.subplots (figsize=figsize)
            _ = ax.plot (x_scale, signal, color=signal_color, linewidth=signal_linewidth, zorder=1)

            legend = []
            #Plot Kmer boundaries if required and available
            if "basecall" in self.data:
                ymin, ymax = ax.get_ylim()
                y1, y2 = ymax, ymax+((ymax-ymin)/10)

                for row in self.data["basecall"].kmers:
                    if row["end"] != row["start"] and row["end"] >= start and row["start"] <= end:
                        _ = ax.hlines (y=row["median"], xmin=row["start"], xmax=row["end"], linewidth=2, color="red", zorder=2, alpha=0.75)
                legend.append (mpatches.Patch(color='red', label='Basecall'))

            # Define title
            title = "Read ID {} / Raw:{:,}".format (self.read_id, len(signal))
            _ = ax.set_title (title)
            _ = ax.set_xlim(start, end)
            if legend:
                _ = ax.legend (handles=legend, frameon=False, ncol=len(legend))
        return fig, ax

    def _to_db (self, grp):
        """Write object into an open h5 group"""

        # Save Metadata
        md_grp = grp.create_group("metadata")
        write_attrs (md_grp.create_group("context_tags"), self.metadata["context_tags"])
        write_attrs (md_grp.create_group("tracking_id"), self.metadata["tracking_id"])
        write_attrs (md_grp.create_group("channel_id"), self.metadata["channel_id"])

        # Save raw
        if "raw" in self.data:
            raw_grp = grp.create_group("raw")
            self.data["raw"]._to_db (raw_grp)

        #Save basecall
        if "basecall" in self.data:
            basecall_grp = grp.create_group("basecall")
            self.data["basecall"]._to_db (basecall_grp)
            # Add hard link to raw in basecall
            basecall_grp["signal"] = raw_grp["signal"]

    #~~~~~~~~~~~~~~CLASS METHODS~~~~~~~~~~~~~~#
    @classmethod
    def from_db (cls, read_id, grp):

        # Get Metadata values
        metadata = OrderedDict ()
        data = OrderedDict ()

        for k, g in grp.items():
            # Get metadata
            if k == "metadata":
                metadata["context_tags"] = parse_attrs (g.get("context_tags"))
                metadata["tracking_id"] = parse_attrs (g.get("tracking_id"))
                metadata["channel_id"] = parse_attrs (g.get("channel_id"))
            elif k == "raw":
                data[k] = Raw.from_db (g)
            elif k == "basecall":
                data[k] = Basecall.from_db (g)
            # elif k.startswith ("alignment"):
            #     data[k] = Alignment.from_db (g)
            # elif k.startswith ("eventalign"):
            #     data[k] = Eventalign.from_db (g)

        return Fast5 (read_id=read_id, metadata=metadata, data=data)

    @classmethod
    def from_fast5 (cls,
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
        data = OrderedDict()
        metadata = OrderedDict ()

        with h5py.File(fast5_fn, "r") as f:

            # Get Metadata values
            try:
                metadata["context_tags"] = parse_attrs (f.get("UniqueGlobalKey/context_tags"))
                metadata["tracking_id"] = parse_attrs (f.get("UniqueGlobalKey/tracking_id"))
                metadata["channel_id"] = parse_attrs (f.get("UniqueGlobalKey/channel_id"))
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Metadata values")

            # Get raw values
            try:
                raw_grp = list(f['/Raw/Reads'].values())[0]
                data ["raw"] = Raw.from_fast5 (grp=raw_grp, signal_normalization=signal_normalization)
                read_id = raw_grp.attrs["read_id"].decode("utf8")
            except (KeyError, IndexError, TypeError) as E:
                raise Fast5Error ("No Raw Value")

            # Get basecall values if available
            try:
                basecall_group = f["Analyses"][basecall_id]
                data ["basecall"] = Basecall.from_fast5 (grp=basecall_group, signal=data["raw"].signal)
            except (KeyError, IndexError, TypeError) as E:
                if basecall_required:
                    raise Fast5Error ("No Basecall Value")

        return Fast5 (read_id=read_id, metadata=metadata, data=data)
