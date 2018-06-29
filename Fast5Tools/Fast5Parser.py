#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import multiprocessing as mp
import sys
import os
from glob import iglob
from time import time
from collections import Counter
import argparse
import shelve

# Third party imports
import numpy as np

# Local imports
from Fast5Tools import __version__
from Fast5Tools.Helper_fun import stdout_print
from Fast5Tools.Fast5 import Fast5, Fast5Error

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Fast5Parser ():
    """
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#

    def __init__ (self,
        fast5_dir,
        output_db_file,
        basecall_group='Basecall_1D_000',
        raw_read_num=0,
        error_on_missing_basecall=True,
        threads = 4,
        verbose = False,
        **kwargs):
        """
        Verify some of the critical parameters and parse the fast5 files using multiple threads.
        Fast5 object are save in a python Shelve for later use
        * fast5_dir: STR
            Path to the folder containing Fast5 files (can be in multiple subfolder)
        * output_db_file: STR
            Path to a directory, where to store the results and logs
        * basecall_group: STR (default 'Basecall_1D_000')
            Name of the analysis group in the fast5 file containing the basecalling information
        * raw_read_num: INT (default 0)
            Index of the read in the Raw group
        * threads: INT (default 4)
            Total number of threads. Minimum = 4
        """
        if verbose: stdout_print ("Initialise Karystos\n")

        # Tests values
        assert os.path.isdir (fast5_dir) and os.access (fast5_dir, os.R_OK), ("Cannot read in the source fast5 folder")
        assert threads >= 4, ("threads should be greater that 4")

        # Fast5 option
        self.fast5_dir = os.path.abspath(fast5_dir)
        self.basecall_group = basecall_group
        self.raw_read_num = raw_read_num

        # Other options
        self.threads = threads
        self.fast5_worker_threads = threads-2
        self.verbose = verbose
        self.error_on_missing_basecall = error_on_missing_basecall
        self.output_db_file = output_db_file

    def __call__ (self):
        """
        All you cores belong to us!
        """
        if self.verbose: stdout_print ("Parse Fast5 files and find blocks\n")

        # Init Multiprocessing variables
        manager = mp.Manager()
        fast5_fn_q = manager.Queue() # Queue for fast5 file names
        fast5_obj_q = manager.Queue() # Queue for blocks found

        # Init processes for file reading, distributed parsing and writing
        ps_list = []
        ps_list.append (mp.Process (target=self._fast5_list_worker, args=(fast5_fn_q,)))
        for i in range (self.fast5_worker_threads):
            ps_list.append (mp.Process (target=self._fast5_parse_worker, args=(fast5_fn_q, fast5_obj_q)))
        ps_list.append (mp.Process (target=self._write_db_worker, args=(fast5_obj_q,)))

        # Start processes and blocks until the process is finished
        for ps in ps_list:
            ps.start()
        for ps in ps_list:
            ps.join()

    def __repr__(self):
        msg = "[{}]\n".format (self.__class__.__name__)

        # list all values in object dict in alphabetical order
        for k, v in sorted (self.__dict__.items(), key=lambda t: t[0]):
            msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _fast5_files_generator (self, fast5_dir, **kwargs):
        """
        Generator returning fast5 files found recursively starting from a given folder.
        The recursivity stops as soon as a file matching the extension is found.
        * fast5_dir: STR
            Path to the folder containing Fast5 files (can be in multiple subfolder)
        """
        # In the case where the folder is a file
        if os.path.isdir(fast5_dir):

            # If matching files in the folder
            fast5_found=False
            for fast5 in iglob (os.path.join(fast5_dir, "*.fast5")):
                yield fast5
                fast5_found=True

            # If no matching file go deeper until a leaf containing fast5 is found
            if not fast5_found:
                for item in os.listdir(fast5_dir):
                    for fast5 in self._fast5_files_generator (os.path.join(fast5_dir, item)):
                        yield fast5

    def _fast5_list_worker (self, fast5_fn_q):
        """
        Mono-threaded worker adding fast5 file found througout a directory tree
        to a feeder queue for the multiprocessing processing workers
        """
        # Load an input Queue with fast5 file path
        for fast5_fn in self._fast5_files_generator (fast5_dir = self.fast5_dir):
            fast5_fn_q.put(fast5_fn)

        # Add 1 poison pill per worker thread
        for i in range (self.fast5_worker_threads):
            fast5_fn_q.put(None)

    def _fast5_parse_worker (self, fast5_fn_q, fast5_obj_q):
        """
        Multi-threaded workers in charge of parsing fast5 file. From valid fast5 files, block matching barcode
        geometries are extracted and add to an output shared memory list. In the same time, metadata are
        collected to be used later to generate metrics
        """
        for fast5_fn in iter(fast5_fn_q.get, None):

            # Try to create a Fast5 object from the Fast5 file path
            try:
                f = Fast5 (
                    fast5_fn = fast5_fn,
                    basecall_group = self.basecall_group,
                    raw_read_num = self.raw_read_num,
                    error_on_missing_basecall = self.error_on_missing_basecall)
                fast5_obj_q.put (f)

            # If an error happened just put it in the out queue
            except Fast5Error as E:
                fast5_obj_q.put (E)

        # Add poison pill in queues
        fast5_obj_q.put (None)

    def _write_db_worker (self, fast5_obj_q):
        """
        Save the block object found in a shelves database while emptying the Queue
        """
        # Init Counters
        err_counter = Counter ()
        n_invalid = n_valid = 0
        t = time()

        # Create shelves database to store the blocks
        with shelve.open (self.output_db_file, flag = "n") as db:
            for _ in range (self.fast5_worker_threads):
                for item in iter (fast5_obj_q.get, None):

                    # If error at instantiation
                    if isinstance (item, Fast5Error):
                        err_counter[item.err_msg] +=1
                        n_invalid +=1

                    # Add new entry in the database
                    elif isinstance (item, Fast5):
                        n_valid += 1
                        db[item.read_id] = item

                    if time()-t >= 0.1:
                        stdout_print("Valid files:{:,} Invalid File:{:,}{}\r".format (n_valid, n_invalid, " "*10))
                        t = time()

        stdout_print("Valid files:{:,} Invalid File:{:,}{}\n".format (n_valid, n_invalid, " "*10))
        if err_counter:
            stdout_print ("Invalid fast5 files summary\n")
            for i, j in err_counter.items():
                stdout_print ("\t{}:{:,}\n".format(i,j))

#~~~~~~~~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~~~~~~~~#

def main ():
    # Define parser options
    parser = argparse.ArgumentParser(description="""Go get those Fast5!""")
    parser.add_argument('--version', '-v', action='version', version=__version__)
    # Required args
    parser.add_argument("-i", "--fast5_dir", required=True, help="Path to the folder containing Fast5 files")
    parser.add_argument("-o", "--output_db_file", required=True, help="Path to the output database file")
    # Optional args
    parser.add_argument("--basecall_group", default='Basecall_1D_000', help="Name of the analysis group in the fast5 file containing the basecalling information")
    parser.add_argument("--raw_read_num", default=0, type=int, help="Index of the read in the Raw")
    parser.add_argument("--threads", default=4, type=int, help="Total number of threads. Minimum = 4")
    parser.add_argument("--verbose", default=False, type=bool, help="If True will be more chatty")
    a = parser.parse_args()

    f = Fast5Parser(
        fast5_dir = a.fast5_dir,
        output_db_file=a.output_db_file,
        basecall_group = a.basecall_group,
        raw_read_num = a.raw_read_num,
        threads = a.threads,
        verbose = a.verbose)
    f()
