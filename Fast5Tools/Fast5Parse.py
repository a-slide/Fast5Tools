  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import multiprocessing as mp
from time import time
from collections import Counter

# Third party imports
import h5py
import numpy as np

# Local imports
from Fast5Tools.Helper_fun import stderr_print, recursive_file_gen, access_dir
from Fast5Tools.Fast5 import Fast5, Fast5Error
from Fast5Tools.Fast5Wrapper import Fast5Wrapper

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
def Fast5Parse (
    fast5_dir,
    db_file,
    basecall_id='Basecall_1D_000',
    signal_normalization="zscore",
    threads = 4,
    max_fast5=None,
    verbose = False,
    **kwargs):
    """
    Verify some of the critical parameters and parse the fast5 files using multiple threads.
    Fast5 object are save in a python Shelve for later use
    * fast5_dir: STR
        Path to the folder containing Fast5 files (can be in multiple subfolder)
    * db_file: STR
        Path to a directory, where to store the results and logs
    * basecall_id: STR (default 'Basecall_1D_000')
        Name of the analysis group in the fast5 file containing the basecalling information
    * signal_normalization (default 'zscore')
        Normalization strategy of the raw signal. Can be None or 'zscore'
    * threads: INT (default 2)
        Total number of threads. Minimum = 2
    """
    if verbose: stderr_print ("Initialise\n")

    # Tests values
    assert access_dir (fast5_dir), "Cannot read in the source fast5 folder"
    assert threads >= 2, "threads should be greater that 2"

    if verbose: stderr_print ("Parse Fast5 files and save to database\n")

    # Init Multiprocessing variables
    manager = mp.Manager()
    fast5_fn_q = manager.Queue(maxsize=1000) # Queue for fast5 file names
    fast5_obj_q = manager.Queue(maxsize=1000) # Queue for blocks found

    # list_worker process
    fl_ps = mp.Process (
        target=_fast5_list_worker,
        args=(fast5_fn_q, fast5_dir, threads, max_fast5))

    # fast5_parse_worker processes
    fp_ps_list = []
    for i in range (threads):
        fp_ps_list.append (mp.Process (
            target=_fast5_parse_worker,
            args=(fast5_fn_q, fast5_obj_q, basecall_id, signal_normalization)))

    # write_db_worker process
    wd_ps = mp.Process (
        target=_write_db_worker,
        args=(fast5_obj_q, db_file, threads, verbose))

    # Start processes
    fl_ps.start ()
    for ps in fp_ps_list:
        ps.start ()
    wd_ps.start ()

    # Join processes
    fl_ps.join ()
    for ps in fp_ps_list:
        ps.join ()
    wd_ps.join ()

    # Return Fast5Wrapper object for further
    #return Fast5Wrapper (db_file, verbose) ###################################################

#~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

def _fast5_list_worker (fast5_fn_q, fast5_dir, threads, max_fast5):
    """
    Mono-threaded worker adding fast5 file found througout a directory tree
    to a feeder queue for the multiprocessing processing workers
    """
    # Load an input Queue with fast5 file path

    for i, fast5_fn in enumerate(recursive_file_gen (dir=fast5_dir, ext="fast5")):
        if max_fast5 and i == max_fast5:
            break
        fast5_fn_q.put(fast5_fn)

    # Add 1 poison pill per worker thread
    for i in range (threads):
        fast5_fn_q.put(None)

def _fast5_parse_worker (fast5_fn_q, fast5_obj_q, basecall_id, signal_normalization):
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
                basecall_id = basecall_id,
                signal_normalization = signal_normalization)
            fast5_obj_q.put (f)

        # If an error happened just put it in the out queue
        except Fast5Error as E:
            fast5_obj_q.put (E)

    # Add poison pill in queues
    fast5_obj_q.put (None)

def _write_db_worker (fast5_obj_q, db_file, threads, verbose):
    """
    Save the block object found in a shelves database while emptying the Queue
    """
    # Init Counters
    err_counter = Counter ()
    n_invalid = n_valid = 0
    t = time()
    read_id_list = []

    # Create shelves database to store the blocks
    with h5py.File(db_file, "w") as fp:

        all_fast5_grp = fp.create_group("fast5")

        for _ in range (threads):
            for item in iter (fast5_obj_q.get, None):

                # If error at instantiation
                if isinstance (item, Fast5Error):
                    err_counter[item.err_msg] +=1
                    n_invalid +=1

                # Add new entry in the database
                elif isinstance (item, Fast5):
                    n_valid += 1
                    read_id = item.read_id
                    fast5_grp = all_fast5_grp.create_group(read_id)
                    item._to_hdf5 (grp = fast5_grp)
                    read_id_list.append (read_id)

                    #### Collect medatata for global metadata

                if time()-t >= 0.2:
                    if verbose: stderr_print("\tValid files:{:,} Invalid File:{:,}\r".format (n_valid, n_invalid))
                    t = time()

        # Write metadata
        all_fast5_grp.attrs.create ("valid_fast5", n_valid)
        all_fast5_grp.attrs.create ("invalid_fast5", n_invalid)

        # Write read_ids:
        read_id_arr = np.array (read_id_list, dtype="<S40")
        all_fast5_grp.create_dataset ("read_ids", data=read_id_arr)

    stderr_print("\tValid files:{:,} Invalid File:{:,}\n".format (n_valid, n_invalid))
    if err_counter:
        stderr_print ("\tInvalid fast5 files summary\n")
        for i, j in err_counter.items():
             stderr_print ("\t\t{}:{:,}\n".format(i,j))
