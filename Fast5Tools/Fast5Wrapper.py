  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
from itertools import islice
from time import time
import sys

# Third party imports
import numpy as np
import pysam
import h5py

# Local imports
from Fast5Tools.Helper_fun import stderr_print, access_file
from Fast5Tools.Fast5 import Fast5, Fast5Error

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Fast5Wrapper ():
    """
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__ (self,
        db_fn,
        verbose = False,
        **kwargs):

        self.db_fn = db_fn
        self.verbose = verbose

        # Check if db is readable
        if not access_file (self.db_fn):
            stderr_print ("Can not access or read the database file\n")
            sys.exit()

        # If index not accessible or not readeable
        if self.verbose:
            stderr_print ("Load reads ids\n")

        with h5py.File(self.db_fn, "r+") as db:
            self.read_id_list = db["read_ids"].value.astype("<U40")

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.db_fn)
        m+= "\tNumber of reads:{:,}".format(len(self))
        return (m)

    def __enter__(self):
        if self.verbose:
            stderr_print ("Open hdf5 database\n")
        self.db = h5py.File(self.db_fn, "r+")
        return self

    def __exit__(self, type, value, traceback):
        if self.verbose:
            stderr_print ("Close hdf5 database\n")
        self.db.close()

        def __len__(self):
            return len(self.read_id_list)

    #~~~~~~~~~~~~~~ GETTER METHODS ~~~~~~~~~~~~~~#

    def get_fast5 (self, read_id):
        grp = self.db["fast5"][read_id]
        return Fast5.from_db (read_id=read_id, grp=grp)

    def sample_fast5 (self, n=10):
        l =[]
        for read_id in np.random.choice (self.read_id_list, n):
            grp = self.db["fast5"][read_id]
            l.append (Fast5.from_db (read_id=read_id, grp=grp))
        return l

    def iter_fast5 (self):
        for read_id in self.read_id_list:
            grp = self.db["fast5"][read_id]
            yield (Fast5.from_db (read_id=read_id, grp=grp))

    #~~~~~~~~~~~~~~MAIN PUBLIC METHODS~~~~~~~~~~~~~~#
    def add_bam_alignment (self,
        fn,
        analysis_name = None,
        only_primary = False,
        verbose=False):
        """
        Parse a bam/sam file
        """
        c = Counter ()
        t = time ()

        if verbose:
            stderr_print ("Parse alignment file {}\n".format(fn))

        # Create root analysis group
        full_analysis_name = "alignment_{}".format(analysis_name) if analysis_name else "alignment"
        assert not full_analysis_name in self.db, "This analysis group already exists. Enter a different suffix"
        analysis_grp = self.db.create_group (full_analysis_name)

        with pysam.AlignmentFile (fn) as bam:

            if verbose:
                stderr_print ("\tParse bam header\n")

            t = time ()
            for r in bam:

                # Counter update
                if time()-t >= 0.2:
                    if verbose:
                        stderr_print ("\tValid hits:{:,}\tPrimary:{:,}\tSecondary:{:,}\tNot in database:{:,}\tUnmapped:{:,}\r".format (
                            c["valid_hit"], c["primary_hit"], c["secondary_hit"], c["not_in_db"], c["unmapped_hit"]))
                    self.db.flush()
                    t = time()

                # Filter out unmapped and reads not in db
                if r.is_unmapped:
                    c["unmapped_hit"] +=1
                elif r.query_name not in self.read_id_list:
                    c["not_in_db"] +=1

                # Valid hits
                else:
                    if r.is_secondary or r.is_supplementary:
                        c["secondary_hit"] +=1
                        if only_primary:
                            continue
                    else:
                        c["primary_hit"] +=1

                    c["valid_hit"] += 1
                    hit_id = "hit_{:010}".format(c["valid_hit"])
                    qname = r.query_name
                    rname = r.reference_name

                    # Get destination fast5 group analysis group and create hardlink
                    fast5_hit_grp = self.db.require_group ("/fast5/{}/{}/{}".format (qname, full_analysis_name, hit_id))
                    ref_hit_group = analysis_grp.require_group (rname)
                    ref_hit_group[hit_id] = fast5_hit_grp

                    # Add hit to alignment
                    fast5_hit_grp.create_dataset ("hit", data=str.encode(r.to_string()))
                    fast5_hit_grp.attrs.create ("read_id", data=str.encode(qname))
                    fast5_hit_grp.attrs.create ("ref_id", data=str.encode(rname))
                    if r.is_secondary or r.is_supplementary:
                        fast5_hit_grp.attrs.create ("hit_type", data=b"secondary")
                    else:
                        fast5_hit_grp.attrs.create ("hit_type", data=b"primary")


            stderr_print ("\tValid hits:{:,}\tPrimary:{:,}\tSecondary:{:,}\tNot in database:{:,}\tUnmapped:{:,}\n".format (
                c["valid_hit"], c["primary_hit"], c["secondary_hit"], c["not_in_db"], c["unmapped_hit"]))
            self.db.flush()

    def add_nanopolish_eventalign (self,
        fn,
        analysis_name = None,
        verbose=False):
        """
        Parse a bam/sam file
        """
        c = Counter ()
        t = time ()

        if verbose:
            stderr_print ("Parse eventalign file {}\n".format(fn))

        # Create root analysis group
        full_analysis_name = "eventalign_{}".format(analysis_name) if analysis_name else "eventalign"
        assert not full_analysis_name in self.db, "This analysis group already exists. Enter a different suffix"
        analysis_grp = self.db.create_group (full_analysis_name)

        with open (fn, "r") as fp:

            # Check header fields
            hs = fp.readline().rstrip().split("\t")
            if not set(["ref_pos", "ref_kmer", "start_idx", "end_idx"]) == set(hs):
                stderr_print ("The input file must be generated by nanopolish eventalign (--print-read-names' and '--signal-index' options)")
                stderr_print ("In addition the files has to be collapsed by kmers using NanopolishComp Eventalign_collapse")
                sys.exit()

            t = time ()
            valid = False
            for line in fp:

                # Split line
                ls = line.rstrip().split("\t")

                # Counter update
                if time()-t >= 0.2:
                    if verbose:
                        stderr_print ("\tValid hits:{:,}\tNot in database:{:,}\r".format (c["valid_hit"], c["not_in_db"]))
                    self.db.flush()
                    t = time()

                # If found read id separator
                if ls[0] == "#":

                    # Write previous kmer if not first
                    if valid:
                        c["valid_hit"] += 1
                        hit_id = "hit_{:010}".format(c["valid_hit"])

                        # Get destination fast5 group analysis group and create hardlink
                        fast5_hit_grp = fast5_grp.require_group ("{}/{}".format (full_analysis_name, hit_id))
                        ref_hit_group = analysis_grp.require_group (rname)
                        ref_hit_group[hit_id] = fast5_hit_grp

                        # Add kmer array
                        kmers = np.array (kmer_list, dtype=[('seq', 'S5'), ('start', np.uint64), ('end', np.uint64), ('ref_pos', np.uint64), ('median', np.float64)])
                        fast5_hit_grp.create_dataset ("kmers", data=kmers, compression="lzf")
                        fast5_hit_grp.attrs.create ("read_id", data=str.encode(qname))
                        fast5_hit_grp.attrs.create ("ref_id", data=str.encode(rname))

                        # Add hard link to raw in basecall
                        fast5_hit_grp["signal"] = fast5_grp["raw/signal"]

                    # Last line exception (flaged with a single #)
                    if len(ls) == 1:
                        break

                    # Start new kmer table
                    qname = ls[1]
                    rname = ls[2]

                    # Save only if qname is in the db
                    if qname in self.read_id_list:
                        valid = True
                        kmer_list = []
                        fast5_grp = self.db.get("fast5/{}".format(qname))
                        signal = fast5_grp.get("raw/signal").value
                    else:
                        valid = False
                        c["not_in_db"] +=1

                # Append line values to the list
                elif valid:
                    seq = ls[1]
                    start = int(ls[2])
                    end = int(ls[3])
                    ref_pos = int(ls[0])
                    sig = signal [start:end]
                    median = np.median(sig) if len(sig) > 0 else np.nan
                    kmer_list.append ((seq, start, end, ref_pos, median))

        stderr_print ("\tValid hits:{:,}\tNot in database:{:,}\n".format (c["valid_hit"], c["not_in_db"]))
        self.db.flush()
