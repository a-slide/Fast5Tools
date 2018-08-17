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
from Fast5Tools.Helper_fun import stderr_print, access_file, parse_attrs, write_attrs
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
        self.verbose = verbose

        # Check if db is readable
        if not access_file (self.db_file):
            stderr_print ("Can not access or read the database file\n")
            sys.exit()

        # If index not accessible or not readeable
        if self.verbose:
            stderr_print ("Load reads ids\n")

        with h5py.File(self.db_file, "r+") as db:
            self.read_id_list = db["read_ids"].value

    def __repr__(self):
        """ Readable description of the object """
        m="[{}] file:{}\n".format(self.__class__.__name__, self.db_file)
        m+= "\tNumber of reads:{:,}".format(len(self))
        return (m)

    def __enter__(self):
        if self.verbose:
            stderr_print ("Open hdf5 database\n")
        self.db = h5py.File(self.db_file, "r+")
        return self

    def __exit__(self, type, value, traceback):
        if self.verbose:
            stderr_print ("Close hdf5 database\n")
        self.db.close()

        def __len__(self):
            return len(self.read_id_list)

    #~~~~~~~~~~~~~~ GETTER METHODS ~~~~~~~~~~~~~~#

    def get_fast5 (self, read_id):
        f_group = self.db["fast5"][read_id]
        return Fast5.from_db (f_group)

    def sample_fast5 (self, n=10):
        l =[]
        for read_id in np.random.choice (self.read_id_list, n):
            f_group = self.db["fast5"][read_id]
            l.append (Fast5.from_db (f_group))
        return l

    def iter_fast5 (self):
        for read_id in self.read_id_list:
            f_group = self.db["fast5"][read_id]
            yield (Fast5.from_db (f_group))

    #~~~~~~~~~~~~~~MAIN PUBLIC METHODS~~~~~~~~~~~~~~#
    def add_bam_alignment (self,
        alignment_fn,
        alignment_num = 0,
        verbose=False):
        """
        Parse a bam/sam file
        """
        c = Counter ()
        t = time ()
        valid_qname = set()
        sync_buffer = 0

        if verbose:
            stderr_print ("Parse alignment file {}\n".format(alignment_fn))

        align_name = "alignment_{0:02d}".format(alignment_num)
        with pysam.AlignmentFile (alignment_fn) as bam:

            # Check header presence
            assert bam.has_index(), "No index/header found in bam file"

            # Create root alignment group
            if align_name in self.db:
                del self.db[align_name]
            align_grp = self.db.create_group (align_name)

            #
            # # Try to parse bam metadata
            # if verbose:
            #     stderr_print ("\tParse bam header\n")
            # try:
            #     write_attrs (grp, bam.header["PG"][0])
            #     write_attrs (grp, bam.header["HD"])
            # except:
            #     pass

            # Create a reference sequence index
            # if verbose:
            #     stderr_print ("\tCreate reference index\n")
            # for ref_dict in bam.header["SQ"]:
            #     ref_name = ref_dict["SN"]
            #     ref_len = ref_dict["LN"]
            #     all_ref_grp = grp.create_group (ref_name)
            #     all_ref_grp.attrs.create (name="ref_len", data=ref_len)


            if verbose:
                stderr_print ("\tParse bam header\n")

            t = time ()
            for r in bam:

                # Counter update
                if time()-t >= 0.2:
                    if verbose:
                        stderr_print ("\tValid hits:{:,}\tPrimary:{:,}\tSecondary:{:,}\tSupplementary:{:,}\tNot in database:{:,}\tUnmapped:{:,}\r".format (
                            c["valid_hit"], c["primary_hit"], c["secondary_hit"], c["supplementary_hit"], c["not_in_db"], c["unmapped_hit"]))
                    self.db.flush()
                    t = time()

                # Filter out unmapped and reads not in db
                if r.is_unmapped:
                    c["unmapped_hit"] +=1
                elif str.encode(r.query_name) not in self.read_id_list:
                    c["not_in_db"] +=1

                # Valid hits
                else:
                    c["valid_hit"] += 1
                    if r.is_secondary:
                        c["secondary_hit"] +=1
                    elif r.is_supplementary:
                        c["supplementary_hit"] +=1
                    else:
                        c["primary_hit"] +=1

                    qname = r.query_name
                    rname = r.reference_name

                    # Get corresponding hit num and fast5 group
                    hit_id = "hit_{0:010d}".format(c["valid_hit"])
                    fast5_grp = self.db.get(f"/fast5/{qname}/")

                    # Remove previous analysis and create new if needed
                    if qname not in valid_qname and align_name in fast5_grp:
                        del fast5_grp [align_name]
                    fast5_hit_group = fast5_grp.require_group (f"{align_name}/{hit_id}")
                    valid_qname.add (qname)

                    # Add Hard link to root alignment group as well
                    ref_group = align_grp.require_group (f"{rname}")
                    ref_group[hit_id] = fast5_hit_group

                    # Add hit to alignment
                    self._add_hit (r, grp=fast5_hit_group)

            stderr_print ("\tValid hits:{:,}\tPrimary:{:,}\tSecondary:{:,}\tSupplementary:{:,}\tNot in database:{:,}\tUnmapped:{:,}\n".format (
                c["valid_hit"], c["primary_hit"], c["secondary_hit"], c["supplementary_hit"], c["not_in_db"], c["unmapped_hit"]))
            self.db.flush()
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

    def _add_hit (self, r, grp):

        grp.attrs.create ("read_name", data=str.encode(r.query_name))
        grp.attrs.create ("read_start", data=int(r.query_alignment_start))
        grp.attrs.create ("read_end", data=int(r.query_alignment_end))
        grp.attrs.create ("ref_name", data=str.encode(r.reference_name))
        grp.attrs.create ("ref_start", data=int(r.reference_start))
        grp.attrs.create ("ref_end_", data=int(r.reference_end))
        grp.attrs.create ("ref_strand", data= b"-" if r.is_reverse else b"+")
        grp.attrs.create ("align_len", data=int(r.query_alignment_length))
        grp.attrs.create ("mapq", data=int(r.mapping_quality))
        grp.attrs.create ("align_score", data=int(r.get_tag("AS")))
        if r.is_secondary:
            grp.attrs.create ("type", data=b"secondary")
        elif r.is_supplementary:
            grp.attrs.create ("type", data=b"supplementary")
        else:
            grp.attrs.create ("type", data=b"primary")

        grp.create_dataset ("cigar", data=str.encode(r.cigarstring))
