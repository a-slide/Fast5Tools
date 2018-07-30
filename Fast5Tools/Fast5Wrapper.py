  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
from collections import Counter
import shelve
from itertools import islice
import random

# Third party imports
import numpy as np
import pysam
from tqdm import tqdm

# Local imports
from Fast5Tools.Helper_fun import stdout_print
from Fast5Tools.Fast5 import Fast5, Fast5Error
from Fast5Tools.Basecall import Basecall
from Fast5Tools.Alignment import Alignment, Read
from Fast5Tools.Nanopolish import Nanopolish


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
        summary_c = Counter()
        read_id_c = Counter ()

        if self.verbose:
            stdout_print ("Parse reads from file {}\n".format(alignment_fn))

        with shelve.open (self.db_file, flag="w", writeback=True) as db, pysam.AlignmentFile(alignment_fn) as bam:
            for r in tqdm(bam, disable= not self.verbose):

                if r.is_unmapped:
                    summary_c["unmapped"] +=1
                    continue

                if r.is_secondary or r.is_supplementary:
                    summary_c["secondary"] +=1
                    if not include_secondary:
                        continue

                # Create new analyses entry in Fast5 if never saw before (overwrite existing)
                qname = r.query_name
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
            stdout_print ("\tReads found")
            for i, j in summary_c.most_common():
                stdout_print ("\t{}: {}".format(i,j))
            stdout_print ("\n\tUnique Fast5 with alignments {}\n".format(len(read_id_c)))

    def add_nanopolish (self, nanopolish_fn, analysis_name="Nanopolish"):
        """
        Parse a nanopolish event align file
        """
        pass

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
