  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
import shelve
from itertools import islice
import random
from time import time

# Third party imports
import numpy as np
import pysam
from tqdm import tqdm

# Local imports
from Fast5Tools.Helper_fun import stderr_print
from Fast5Tools.Fast5 import Fast5, Fast5Error
from Fast5Tools.Basecall import Basecall
from Fast5Tools.Alignment import Alignment, Read
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
        read_id_c = Counter ()
        summary_c = Counter()
        t = time ()

        if self.verbose:
            stderr_print ("Parse alignment file {}\n".format(alignment_fn))

        with shelve.open (self.db_file, flag="w", writeback=True) as db, pysam.AlignmentFile(alignment_fn) as fp:
            for r in fp:

                # Counter update
                if self.verbose and time()-t >= 0.2:
                    stderr_print("\tValid hits:{:,}\tInvalid hits:{:,}\tSecondary hits:{:,}\tUnmapped reads:{:,}\r".format (
                        summary_c["valid"], summary_c["invalid"], summary_c["secondary"], summary_c["unmapped"]))
                    t = time()

                qname = r.query_name
                if qname not in db:
                    summary_c["not_in_db"] +=1 ################################# COUNT EVERY LINE NOT IN DB not just READS
                    continue

                if r.is_unmapped:
                    summary_c["unmapped"] +=1
                    continue

                if r.is_secondary or r.is_supplementary:
                    summary_c["secondary"] +=1
                    if not include_secondary:
                        continue

                # Create new analyses entry in Fast5 if never saw before (overwrite existing)
                if not qname in read_id_c:
                    db[qname].analyses[analysis_name] = Alignment ()
                read_id_c [qname] +=1

                # Add read to alignment analysis
                try:
                    db[qname].analyses[analysis_name].add_read (
                        Read (
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
                    if self.verbose: summary_c["valid"] +=1

                except TypeError as E:
                    if self.verbose: summary_c["invalid"] +=1

        if self.verbose:
            stderr_print("\tValid hits:{:,}\tInvalid hits:{:,}\tSecondary hits:{:,}\tUnmapped reads:{:,}\tReads not in db:{:,}\n".format (
                summary_c["valid"], summary_c["invalid"], summary_c["secondary"], summary_c["unmapped"], summary_c["not_in_db"]))
            stderr_print ("\tUnique Fast5 with alignments {}\n".format(len(read_id_c)))

    def add_nanopolish_eventalign (self, eventalign_fn, analysis_name="Nanopolish_eventalign"):
        """
        Parse a nanopolish event align file
        """
        array_dtypes =[
            ('seq', '<U5'),
            ('start', np.uint32),
            ('end', np.uint32),
            ('ref_pos', np.uint32),
            ('mean', np.longdouble),
            ('median', np.longdouble),
            ('std', np.longdouble)]

        read_id_c = Counter ()
        summary_c = Counter()
        t = time ()

        if self.verbose:
            stderr_print ("Parse Nanopolish eventalign file {}\n".format(eventalign_fn))

        with shelve.open (self.db_file, flag="w", writeback=True) as db, open (eventalign_fn, "r") as fp:

            # Flush header line
            header = next(fp)

            first = True
            for line in fp:

                # Counter update
                if self.verbose and time()-t >= 0.2:
                    stderr_print ("\tValid reads:{:,}\tReads not in db:{:,}\r".format (
                        summary_c["reads"], summary_c["not_in_db"]))
                    t = time()

                # Extract important fields from the file
                ls = line.rstrip().split("\t")
                if len(ls) != 15:
                    summary_c["invalid_lines"] += 1
                cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname = ls[3] ,ls[2] ,int(ls[13]) ,int(ls[14]) ,int(ls[1]) ,ls[0]

                if cur_qname not in db:
                    summary_c["not_in_db"] +=1

                # First line exception
                elif first:
                    first = False
                    qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname
                    raw = db[qname].get_raw()
                    kmer_list = []

                # Fill in the list for the curent read name
                elif cur_qname == qname:
                    if cur_ref_pos == ref_pos:
                        start = cur_start

                    else:
                        rs = raw [start:end]
                        kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                        seq, start, end, ref_pos = cur_seq, cur_start, cur_end, cur_ref_pos

                else:
                    rs = raw [start:end]
                    kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                    db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
                    qname, seq, start, end, ref_pos, rname = cur_qname, cur_seq, cur_start, cur_end, cur_ref_pos, cur_rname

                    # Start a new kmer list
                    raw = db[qname].get_raw()
                    kmer_list = []
                    summary_c["reads"] += 1

            # Last line exception
            if kmer_list:
                rs = raw [start:end]
                kmer_list.insert (0, (seq, start, end, ref_pos, np.mean(rs), np.median(rs), np.std(rs)))
                db[qname].analyses[analysis_name] = Eventalign (ref_name = rname, kmers = np.array (kmer_list, dtype=array_dtypes))
                summary_c["kmers"] += 1
                summary_c["reads"] += 1

        if self.verbose:
            stderr_print ("\tValid reads:{:,}\tReads not in db:{:,}\r".format (
                summary_c["reads"], summary_c["not_in_db"]))


    #~~~~~~~~~~~~~~PROPERTY HELPER AND MAGIC METHODS~~~~~~~~~~~~~~#
    def head (self, n=5):
        l =[]
        with shelve.open (self.db_file, flag = "r") as db:
            for f in islice( db.values (), n):
                l.append (f)
        return l

    def sample (self, n=10):
        with shelve.open (self.db_file, flag = "r") as db:
            l = []
            for k in random.sample (list(db.keys()), n):
                l.append (db[k])
        return l

    def __len__ (self):
        with shelve.open (self.db_file, flag = "r") as db:
            return (len(db))

    def __iter__ (self):
        with shelve.open (self.db_file, flag = "r") as db:
            for f in db.values ():
                yield (f)

    def __getitem__(self, items):
        with shelve.open (self.db_file, flag = "r") as db:
            if items in db:
                return db[items]
            else:
                return None
