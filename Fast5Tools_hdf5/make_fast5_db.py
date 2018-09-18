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
from Fast5Tools_hdf5.Helper_fun import stderr_print, recursive_file_gen, access_dir
from Fast5Tools_hdf5.Fast5 import Fast5, Fast5Error

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
def make_fast5_db (
    fast5_dir,
    db_fn,
    basecall_id='Basecall_1D_000',
    signal_normalization="zscore",
    threads = 2,
    max_fast5=0,
    basecall_required=False,
    verbose = False,
    **kwargs):
    """
    fast5_parse read all fast5 files in a directory recursively, extract raw signal, metadata and Albacore basecalling values (if available).
    The fast5 objects generated are stored in a HDF5 database
    * fast5_dir: STR
        Path to the folder containing Fast5 files (can be in multiple subfolder)
    * db_fn: STR
        Path to a directory, where to store the results and logs
    * basecall_id: STR (default 'Basecall_1D_000')
        Name of the analysis group in the fast5 file containing the basecalling information
    * signal_normalization (default 'zscore')
        Normalization strategy of the raw signal. Can be None or 'zscore'
    * threads: INT (default 2)
        Total number of threads. Minimum = 2
    * max_fast5 (default = 0)
        Maximum number of file to try to parse. 0 to deactivate
    """
    if verbose: stderr_print ("Initialise\n")

    # Tests values
    assert access_dir (fast5_dir), "Cannot read in the source fast5 folder"
    ############################################################################ Test if dest is writeable
    assert threads >= 2, "threads should be greater that 2"

    if verbose: stderr_print ("Parse Fast5 files and save to database\n")

    # Init Multiprocessing variables
    manager = mp.Manager()
    fast5_fn_q = manager.Queue(maxsize=1000) # Queue for fast5 file names
    fast5_obj_q = manager.Queue(maxsize=1000) # Queue for blocks found

    # list_worker process
    ps_list = []
    ps_list.append (mp.Process ( target=_fast5_list_worker,args=(fast5_fn_q, fast5_dir, threads, max_fast5)))
    for i in range (threads):
        ps_list.append (mp.Process (target=_fast5_parse_worker, args=(fast5_fn_q, fast5_obj_q, basecall_id, signal_normalization, basecall_required)))
    ps_list.append (mp.Process (target=_write_db_worker, args=(fast5_obj_q, db_fn, threads, verbose)))

    try:
        # Start processes
        for ps in ps_list:
            ps.start ()
        # Join processes
        for ps in ps_list:
            ps.join ()

    # Kill processes if early stop
    except (BrokenPipeError, KeyboardInterrupt) as E:
        if verbose: stderr_print ("Early stop. Kill processes\n")
        for ps in ps_list:
            ps.terminate ()

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

def _fast5_parse_worker (fast5_fn_q, fast5_obj_q, basecall_id, signal_normalization, basecall_required):
    """
    Multi-threaded workers in charge of parsing fast5 file. From valid fast5 files, block matching barcode
    geometries are extracted and add to an output shared memory list. In the same time, metadata are
    collected to be used later to generate metrics
    """
    for fast5_fn in iter(fast5_fn_q.get, None):

        # Try to create a Fast5 object from the Fast5 file path
        try:
            f = Fast5.from_fast5 (
                fast5_fn = fast5_fn,
                basecall_id = basecall_id,
                signal_normalization = signal_normalization,
                basecall_required=basecall_required)
            fast5_obj_q.put (f)

        # If an error happened just put it in the out queue
        except Fast5Error as E:
            fast5_obj_q.put (E)

    # Add poison pill in queues
    fast5_obj_q.put (None)

def _write_db_worker (fast5_obj_q, db_fn, threads, verbose):
    """
    Save the block object found in a shelves database while emptying the Queue
    """
    # Init Counters
    err_counter = Counter ()
    n_invalid = n_valid = 0
    t = time()
    read_id_list = []
    buffer = 0

    # Create a hdf5 database to store the Fast5
    with h5py.File(db_fn, "w") as fp:
        all_fast5_grp = fp.create_group("fast5")

        t = time()
        for _ in range (threads):
            for item in iter (fast5_obj_q.get, None):

                # If error at instantiation
                if isinstance (item, Fast5Error):
                    err_counter[item.err_msg] +=1
                    n_invalid +=1

                # Add new entry in the database
                elif isinstance (item, Fast5):
                    n_valid += 1
                    buffer  += 1
                    read_id = item.read_id
                    fast5_grp = all_fast5_grp.create_group(read_id)
                    item._to_db (grp = fast5_grp)
                    read_id_list.append (read_id)

                if time()-t >= 0.2:
                    if verbose:
                        stderr_print("\tValid files:{:,} Invalid File:{:,}\r".format (n_valid, n_invalid))
                    fp.flush()
                    t = time()

        # Write metadata
        all_fast5_grp.attrs.create ("valid_fast5", n_valid)
        all_fast5_grp.attrs.create ("invalid_fast5", n_invalid)

        # Write read_ids:
        read_id_arr = np.array (read_id_list, dtype="<S40")
        fp.create_dataset ("read_ids", data=read_id_arr, compression="lzf")

    stderr_print("\tValid files:{:,} Invalid File:{:,}\n".format (n_valid, n_invalid))
    if err_counter:
        stderr_print ("\tInvalid fast5 files summary\n")
        for i, j in err_counter.items():
             stderr_print ("\t\t{}:{:,}\n".format(i,j))
