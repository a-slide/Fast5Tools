# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np
import pandas as pd

# Local import
from Fast5Tools.Helper_fun import write_attrs

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Raw (object):
    """
    Represent and summarize basecalling informations from Albacore
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, signal, metadata, normalization, **kwargs):
        """
        """
        # Self variables
        self.signal = signal
        self.metadata = metadata

        if normalization == "zscore":
            self.signal = (self.signal - self.signal.mean()) / self.signal.std()
            self.metadata ["normalization"] = "zscore"

    def __repr__(self):
        """ Readable description of the object """
        m = "[{}]  Signal: {}... / Length: {}".format (self.__class__.__name__, self.signal[0:5], len(self))
        if "normalization" in self.metadata:
            m +=" / Normalization: {}".format(self.metadata ["normalization"])
        return m

    def __len__ (self):
        return len(self.signal)

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#
    @property
    def to_series (self):
        return pd.Series (self.signal)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def get_signal (self, start=None, end=None, smoothing_win_size=0):
        """
        * start INT
            If defined the data will start at that value
        * end INT
            If defined the data will end at that value
        * smoothing_win_size INT
            If larger than 0 will smooth the signal with a moving median window of size X
        """

        if smoothing_win_size:
            signal = self._signal_smoothing (win_size=smoothing_win_size)
            return signal [start:end]

        else:
            return self.signal [start:end]

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _signal_smoothing (self, win_size=3, **kwargs):
        """ Smooth the signal using a moving median window.
        """
        # Create an empty array
        signal = np.empty (dtype=self.signal.dtype, shape=self.signal.shape)
        # Iterate window by window over the signal value array and compute the median for each
        for i, j in enumerate (np.arange (0, len(self))):
            signal [i] = np.median (self.signal[j:j+win_size])
        return signal

    def _to_hdf5 (self, grp):
        """Write object into an open h5 group"""
        # Save Metadata
        write_attrs (grp, self.metadata)
        # Save Signal
        grp.create_dataset("signal", data=self.signal)
