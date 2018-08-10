  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
import shelve
from itertools import islice
import random
from time import time
import pickle
import sys

# Third party imports
import numpy as np
import pysam

# Local imports
from Fast5Tools.Helper_fun import stderr_print, access_file
from Fast5Tools.Fast5 import Fast5, Fast5Error
from Fast5Tools.Basecall import Basecall
from Fast5Tools.Alignment import Alignment, Hit
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

        self.verbose = verbose
        self.db_file = db_file
        self.db_index = db_file+".dbi"

        # Check if db is readable
        if not access_file (self.db_file):
            stderr_print ("Can not access or read the database file")
            sys.exit()

        # If index not accessible or not readeable
        if not access_file (self.db_index):
            if self.verbose:
                stderr_print ("Can not access or read the database index. Indexing database")
            # List all read_id keys in db
            with shelve.open (self.db_file, flag = "r") as db:
                self.read_id_list = list(db.keys())
            # Write in a pickled file for next time
            with open (self.db_index, "wb") as fh:
                pickle.dump (self.read_id_list, fh)

        # If index is readeable. Unpickle the read_id list
        else:
            if self.verbose:
                stderr_print ("Load database index")
            with open (self.db_index, "rb") as fh:
                self.read_id_list = pickle.load(fh)

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.db_file)
        m+= "\tNumber of sequences:{:,}".format(len(self))
        return (m)

    #~~~~~~~~~~~~~~MAIN PUBLIC METHODS~~~~~~~~~~~~~~#
    def add_bam_alignment (self,
        alignment_fn,
        analysis_name = "Alignment",
        include_secondary = False):
        """
        Parse a bam/sam file
        """
        not_in_db_read_id = set()
        valid_read_id = set()
        valid_hits = invalid_hits = 0
        t = time ()

        if self.verbose:
            stderr_print ("Parse alignment file {}\n".format(alignment_fn))

        with shelve.open (self.db_file, flag="w") as db, pysam.AlignmentFile(alignment_fn) as fp:
            for r in fp:

                # Counter update
                if self.verbose and time()-t >= 0.2:
                    stderr_print("\tValid reads:{:,}\tValid hits:{:,}\tReads not in database:{:,}\tSkiped unmapped and secondary:{:,}\r".format (
                        len(valid_read_id), valid_hits, len(not_in_db_read_id), invalid_hits))
                    t = time()

                qname = r.query_name

                if qname not in db:
                    not_in_db_read_id.add(qname)

                elif r.is_unmapped or (not include_secondary and (r.is_secondary or r.is_supplementary)):
                    invalid_hits +=1

                else:
                    # Create new analyses entry in Fast5 if never saw before (overwrite existing)
                    if not qname in valid_read_id:
                        db[qname].analyses[analysis_name] = Alignment ()

                    # Add read to alignment analysis
                    db[qname].analyses[analysis_name].add_read (
                        Hit (
                            qname = qname,
                            qstart = int (r.query_alignment_start),
                            qend = int (r.query_alignment_end),
                            qlen = int (r.query_length),
                            rname = r.reference_name,
                            rstart = int (r.reference_start),
                            rend = int (r.reference_end),
                            rlen = int (r.reference_length),
                            strand = "-" if r.is_reverse else "+",
                            align_len = int (r.query_alignment_length),
                            mapq = int (r.mapping_quality),
                            align_score = int (r.get_tag("AS"))))
                    valid_read_id.add(qname)
                    valid_hits +=1

        if self.verbose:
            stderr_print("\tValid reads:{:,}\tValid hits:{:,}\tReads not in database:{:,}\tSkiped unmapped and secondary:{:,}\n".format (
                len(valid_read_id), valid_hits, len(not_in_db_read_id), invalid_hits))

    def add_nanopolish_eventalign (self,
        eventalign_fn,
        analysis_name="Nanopolish_eventalign"):
        """
        Parse a nanopolish event align file
        """
        array_dtypes =[
            ('seq', '<U5'),
            ('start', np.uint32),
            ('end', np.uint32),
            ('ref_pos', np.uint32),
            ('mean', np.float64),
            ('median', np.float64),
            ('std', np.float64)]

        not_in_db_read_id = set()
        valid_read_id = set()
        valid_kmers = empty_kmers = 0
        t = time ()

        if self.verbose:
            stderr_print ("Parse Nanopolish eventalign file {}\n".format(eventalign_fn))

        with shelve.open (self.db_file, flag="w") as db, open (eventalign_fn, "r") as fp:

            # Check header fields
            header = next(fp)
            hs = header.rstrip().split("\t")
            if not "read_name" in hs:
                raise ValueError ("Nanopolish eventalign has to be ran with '--print-read-names' option")
            if not "start_idx" in hs or not "end_idx" in hs:
                raise ValueError ("Nanopolish eventalign has to be ran with '--signal-index' option")
            if len(hs) != 15:
                raise ValueError ("Missing fields in the nanopolish eventalign")

            first = True
            for line in fp:

                # Counter update
                if self.verbose and time()-t >= 0.2:
                    stderr_print("\tValid reads:{:,}\tValid kmers:{:,}\tEmpty_kmers:{:,}\tReads not in database:{:,}\r".format (
                        len(valid_read_id), valid_kmers, empty_kmers, len(not_in_db_read_id)))
                    t = time()

                # Extract important fields from the file
                ls = line.rstrip().split("\t")
                cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname = ls[3] ,ls[2] ,int(ls[13]) ,int(ls[14]) ,int(ls[1]) ,ls[0]

                if cur_qname not in db:
                    not_in_db_read_id.add(cur_qname)

                # First line exception
                elif first:
                    first = False
                    # Start a new kmer list
                    qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname
                    raw = db[qname].get_raw()
                    kmer_list = []

                # Fill in the list for the curent read name
                elif cur_qname == qname:
                    if cur_ref_pos == ref_pos:
                        # Update start position (Nanopolish coord are reverse compared with raw signal)
                        start = cur_start

                    else:
                        # Add new kmer
                        rs = raw [start:end]
                        if rs.size ==0 :
                            kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
                            empty_kmers +=1
                        else:
                            kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                            valid_kmers +=1
                        # Start new kmer
                        seq, start, end, ref_pos = cur_seq, cur_start, cur_end, cur_ref_pos
                        # Update counters
                        valid_kmers +=1

                else:
                    # Add new kmer and add read_name analysis entry
                    rs = raw [start:end]
                    if rs.size ==0 :
                        kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
                        empty_kmers +=1
                    else:
                        kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                        valid_kmers +=1
                    db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
                    # Start a new kmer list
                    qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname
                    raw = db[qname].get_raw()
                    kmer_list = []
                    # Update counters
                    valid_read_id.add(qname)
                    valid_kmers +=1

            # Last line exception
            if kmer_list:
                # Add new kmer and add read_name analysis entry
                rs = raw [start:end]
                if rs.size ==0 :
                    kmer_list.insert (0, (seq, start, end, ref_pos, np.nan, np.nan, np.nan))
                    empty_kmers +=1
                else:
                    kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                    valid_kmers +=1
                db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
                # Update counters
                valid_read_id.add(qname)

        if self.verbose:
            stderr_print("\tValid reads:{:,}\tValid kmers:{:,}\tEmpty_kmers:{:,}\tReads not in database:{:,}\n".format (
                len(valid_read_id), valid_kmers, empty_kmers, len(not_in_db_read_id)))

    #~~~~~~~~~~~~~~PROPERTY HELPER AND MAGIC METHODS~~~~~~~~~~~~~~#
    def head (self, n=5):
        l =[]
        with shelve.open (self.db_file, flag = "r") as db:
            for read_id in self.read_id_list [:n]:
                f = db[read_id]
                l.append (f)
        return l

    def sample (self, n=10):
        l =[]
        with shelve.open (self.db_file, flag = "r") as db:
            for read_id in random.sample (self.read_id_list, n):
                f = db[read_id]
                l.append (f)
        return l

    def __len__ (self):
        return (len(self.read_id_list))

    def __iter__ (self):
        with shelve.open (self.db_file, flag = "r") as db:
            for read_id in self.read_id_list:
                f = db[read_id]
                yield (f)

    def __getitem__(self, items):
        with shelve.open (self.db_file, flag = "r") as db:
            if items in self.read_id_list:
                return db[items]
            else:
                return None

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
