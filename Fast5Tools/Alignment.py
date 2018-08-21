# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
import shelve

# Third party imports
import pysam

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Hit (object):
    """
    """
    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __init__(self, r, **kwargs):
        """Extract fields from pysam AlignedSegment object they are not pickeable"""
        self.qname = r.query_name
        self.qstart = int (r.query_alignment_start)
        self.qend = int (r.query_alignment_end)
        self.qlen = int (r.query_length)
        self.rname = r.reference_name
        self.rstart = int (r.reference_start)
        self.rend = int (r.reference_end)
        self.rlen = int (r.reference_length)
        self.strand = "-" if r.is_reverse else "+"
        self.align_len = int (r.query_alignment_length)
        self.mapq = int (r.mapping_quality)
        self.align_score = int (r.get_tag("AS"))

    def __repr__(self):
        """ Readable description of the object """
        return "Query:{}-{}:{} ({} pb) / Reference:{}-{}:{}({}) ({} pb) / Alignment len:{} / Mapq:{} / Align Score:{}".format(
            self.qname,
            self.qstart,
            self.qend,
            self.qlen,
            self.rname,
            self.rstart,
            self.rend,
            self.strand,
            self.rlen,
            self.align_len,
            self.mapq,
            self.align_score)


    # # #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    # def _add_hit (self, r, grp):
    #
    #     grp.attrs.create ("read_name", data=str.encode(r.query_name))
    #     grp.attrs.create ("read_start", data=int(r.query_alignment_start))
    #     grp.attrs.create ("read_end", data=int(r.query_alignment_end))
    #     grp.attrs.create ("ref_name", data=str.encode(r.reference_name))
    #     grp.attrs.create ("ref_start", data=int(r.reference_start))
    #     grp.attrs.create ("ref_end_", data=int(r.reference_end))
    #     grp.attrs.create ("ref_strand", data= b"-" if r.is_reverse else b"+")
    #     grp.attrs.create ("align_len", data=int(r.query_alignment_length))
    #     grp.attrs.create ("mapq", data=int(r.mapping_quality))
    #     grp.attrs.create ("align_score", data=int(r.get_tag("AS")))
    #     if r.is_secondary:
    #         grp.attrs.create ("type", data=b"secondary")
    #     elif r.is_supplementary:
    #         grp.attrs.create ("type", data=b"supplementary")
    #     else:
    #         grp.attrs.create ("type", data=b"primary")
    #
    #     grp.create_dataset ("cigar", data=str.encode(r.cigarstring))
