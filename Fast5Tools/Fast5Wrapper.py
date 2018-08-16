  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
from itertools import islice
import random
from time import time
import sys

# Third party imports
import numpy as np
import pysam
import h5py

# Local imports
from Fast5Tools.Helper_fun import stderr_print, access_file
from Fast5Tools.Fast5 import Fast5, Fast5Error
from Fast5Tools.Basecall import Basecall
from Fast5Tools.Alignment import Alignment
from Fast5Tools.Eventalign import Eventalign

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Fast5Wrapper ():
    """
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__ (self,
        db_file,
        verbose = False,
        **kwargs):

        self.db_file = db_file

        # Check if db is readable
        if not access_file (self.db_file):
            stderr_print ("Can not access or read the database file\n")
            sys.exit()

        # If index not accessible or not readeable
        with h5py.File(self.db_file, "r") as db:
            self.read_id_list = db["read_ids"].value

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.db_file)
        m+= "\tNumber of reads:{:,}".format(len(self))
        return (m)

    #~~~~~~~~~~~~~~PROPERTY HELPER AND MAGIC METHODS~~~~~~~~~~~~~~#
    def __len__(self):
        return len(self.read_id_list)

    # def get_fast5 (self, read_id, proxy=False):
    #     wit


    # def head (self, n=5):
    #     l =[]
    #     with h5py.File("pass.h5", "r") as db:
    #         fast5_grp = db["fast5"]
    #         for read_id in self.read_id_list [:n]:
    #             f = Fast5Proxy (fast5_grp[read_id])
    #             l.append (f)
    #     return l
    #
    # def sample (self, n=10):
    #     l =[]
    #     with shelve.open (self.db_file, flag = "r") as db:
    #         for read_id in random.sample (self.read_id_list, n):
    #             f = db[read_id]
    #             l.append (f)
    #     return l
    #
    # def __iter__ (self):
    #     with shelve.open (self.db_file, flag = "r") as db:
    #         for read_id in self.read_id_list:
    #             f = db[read_id]
    #             yield (f)
    #
    # def __getitem__(self, items):
    #     with shelve.open (self.db_file, flag = "r") as db:
    #         if items in self.read_id_list:
    #             return db[items]
    #         else:
    #             return None

    #~~~~~~~~~~~~~~MAIN PUBLIC METHODS~~~~~~~~~~~~~~#
    # def add_bam_alignment (self,
    #     alignment_fn,
    #     analysis_name = "Alignment",
    #     max_sync_buffer = 100,
    #     verbose=False):
    #     """
    #     Parse a bam/sam file
    #     """
    #     not_in_db_read_id = set()
    #     valid_read_id = set()
    #     unmapped_read = secondary_hit = 0
    #     t = time ()
    #     sync_buffer = 0
    #
    #     if verbose:
    #         stderr_print ("Parse alignment file {}\n".format(alignment_fn))
    #
    #     with shelve.open (self.db_file, flag="w", writeback=True) as db, pysam.AlignmentFile(alignment_fn) as fp:
    #         for r in fp:
    #
    #             # Counter update
    #             if verbose and time()-t >= 0.2:
    #                 stderr_print("\tValid reads:{:,}\tReads not in database:{:,}\tReads unmapped:{:,}\tSecondary hits:{:,}\r".format (
    #                     len(valid_read_id), len(not_in_db_read_id), unmapped_read, secondary_hit))
    #                 t = time()
    #
    #             qname = r.query_name
    #
    #             if r.is_secondary or r.is_supplementary:
    #                 secondary_hit +=1
    #             elif r.is_unmapped:
    #                 unmapped_read +=1
    #             elif qname not in self.read_id_list:
    #                 not_in_db_read_id.add(qname)
    #             else:
    #                 db[qname].analyses[analysis_name] = Alignment(r)
    #                 valid_read_id.add(qname)
    #
    #                 sync_buffer += 1
    #                 if sync_buffer == max_sync_buffer:
    #                     db.sync()
    #                     sync_buffer = 0
    #
    #         db.sync()
    #         stderr_print("\tValid reads:{:,}\tReads not in database:{:,}\tReads unmapped:{:,}\tSecondary hits:{:,}\n".format (
    #                 len(valid_read_id), len(not_in_db_read_id), unmapped_read, secondary_hit))
    #
    # def add_nanopolish_eventalign (self,
    #     eventalign_fn,
    #     analysis_name="Nanopolish_eventalign",
    #     verbose=False):
    #     """
    #     Parse a nanopolish event align file
    #     """
    #     array_dtypes =[
    #         ('seq', '<U5'),
    #         ('start', np.uint32),
    #         ('end', np.uint32),
    #         ('ref_pos', np.uint32),
    #         ('mean', np.float64),
    #         ('median', np.float64),
    #         ('std', np.float64)]
    #
    #     not_in_db_read_id = set()
    #     valid_read_id = set()
    #     valid_kmers = empty_kmers = 0
    #     t = time ()
    #
    #     if verbose:
    #         stderr_print ("Parse Nanopolish eventalign file {}\n".format(eventalign_fn))
    #
    #     with shelve.open (self.db_file, flag="w") as db, open (eventalign_fn, "r") as fp:
    #
    #         # Check header fields
    #         header = next(fp)
    #         hs = header.rstrip().split("\t")
    #         if not "read_name" in hs:
    #             raise ValueError ("Nanopolish eventalign has to be ran with '--print-read-names' option")
    #         if not "start_idx" in hs or not "end_idx" in hs:
    #             raise ValueError ("Nanopolish eventalign has to be ran with '--signal-index' option")
    #         if len(hs) != 15:
    #             raise ValueError ("Missing fields in the nanopolish eventalign")
    #
    #         first = True
    #         for line in fp:
    #
    #             # Counter update
    #             if verbose and time()-t >= 0.2:
    #                 stderr_print("\tValid reads:{:,}\tValid kmers:{:,}\tEmpty_kmers:{:,}\tReads not in database:{:,}\r".format (
    #                     len(valid_read_id), valid_kmers, empty_kmers, len(not_in_db_read_id)))
    #                 t = time()
    #
    #             # Extract important fields from the file
    #             ls = line.rstrip().split("\t")
    #             cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname = ls[3] ,ls[2] ,int(ls[13]) ,int(ls[14]) ,int(ls[1]) ,ls[0]
    #
    #             if cur_qname not in db:
    #                 not_in_db_read_id.add(cur_qname)
    #
    #             # First line exception
    #             elif first:
    #                 first = False
    #                 # Start a new kmer list
    #                 qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname
    #                 raw = db[qname].get_raw()
    #                 kmer_list = []
    #
    #             # Fill in the list for the curent read name
    #             elif cur_qname == qname:
    #                 if cur_ref_pos == ref_pos:
    #                     # Update start position (Nanopolish coord are reverse compared with raw signal)
    #                     start = cur_start
    #
    #                 else:
    #                     # Add new kmer
    #                     rs = raw [start:end]
    #                     if rs.size ==0 :
    #                         kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
    #                         empty_kmers +=1
    #                     else:
    #                         kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
    #                         valid_kmers +=1
    #                     # Start new kmer
    #                     seq, start, end, ref_pos = cur_seq, cur_start, cur_end, cur_ref_pos
    #                     # Update counters
    #                     valid_kmers +=1
    #
    #             else:
    #                 # Add new kmer and add read_name analysis entry
    #                 rs = raw [start:end]
    #                 if rs.size ==0 :
    #                     kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
    #                     empty_kmers +=1
    #                 else:
    #                     kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
    #                     valid_kmers +=1
    #                 db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
    #                 # Start a new kmer list
    #                 qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname
    #                 raw = db[qname].get_raw()
    #                 kmer_list = []
    #                 # Update counters
    #                 valid_read_id.add(qname)
    #                 valid_kmers +=1
    #
    #         # Last line exception
    #         if kmer_list:
    #             # Add new kmer and add read_name analysis entry
    #             rs = raw [start:end]
    #             if rs.size ==0 :
    #                 kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
    #                 empty_kmers +=1
    #             else:
    #                 kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
    #                 valid_kmers +=1
    #             db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
    #             # Update counters
    #             valid_read_id.add(qname)
    #
    #     stderr_print("\tValid reads:{:,}\tValid kmers:{:,}\tEmpty_kmers:{:,}\tReads not in database:{:,}\n".format (
    #         len(valid_read_id), valid_kmers, empty_kmers, len(not_in_db_read_id)))



    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
