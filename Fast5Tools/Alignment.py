# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
from collections import Counter
import shelve

# Third party imports
import pysam

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Alignment (object):
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
