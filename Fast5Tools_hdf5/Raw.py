# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports

# Third party imports
import numpy as np

# Local import
from Fast5Tools_hdf5.Helper_fun import write_attrs, parse_attrs

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Raw (object):
    """
    Represent and summarize raw data information
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, signal, metadata, **kwargs):
        """
        """
        # Self variables
        self.signal = signal
        self.metadata = metadata

    def __repr__(self):
        """ Readable description of the object """
        m = "Signal: {}... / Length: {}".format (self.signal[0:5], len(self))
        if "normalization" in self.metadata:
            m +=" / Normalization: {}".format(self.metadata ["normalization"])
        return m

    def __len__ (self):
        return len(self.signal)

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

    def _to_db (self, grp):
        """Write object into an open hdf5 group"""
        # Save Metadata
        write_attrs (grp, self.metadata)
        # Save Signal
        grp.create_dataset("signal", data=self.signal, compression="lzf")

    #~~~~~~~~~~~~~~CLASS METHODS~~~~~~~~~~~~~~#
    @classmethod
    def from_fast5 (cls, grp, signal_normalization):

        # Extract metadata
        metadata = parse_attrs (grp)
        # Extract signal
        signal = grp['Signal'].value
        # Normalise signal if required
        if signal_normalization == "zscore":
            signal = (signal - signal.mean()) / signal.std()
            metadata ["normalization"] = "zscore"

        return Raw (signal=signal, metadata=metadata)

    @classmethod
    def from_db (cls, grp):
        return Raw (
            signal = grp.get("signal").value,
            metadata = parse_attrs (grp))