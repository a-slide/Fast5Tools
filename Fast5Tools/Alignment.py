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
    def __init__(self):
        """"""
        self.read_list = []
        self.nread = 0

    def add_read (self, read, **kwargs):
        """"""
        self.read_list.append (read)
        self.nread += 1

    def best_read (self, based_on = "align_score"):
        """"""
        if self.nread == 0:
            return None

        if self.nread == 1:
            return self.read_list [0]

        if based_on == "align_score":
            return sorted(self.read_list, key=lambda x: x.align_score, reverse=True)[0]

        elif based_on == "mapq":
            return sorted(self.read_list, key=lambda x: x.mapq, reverse=True)[0]

        elif based_on == "align_len":
            return sorted(self.read_list, key=lambda x: x.align_len, reverse=True)[0]

    def __repr__(self):
        """ Readable description of the object """
        m = ""
        for r in self.read_list:
            m +="\t\t{}\n".format(r)
        return (m)

class Hit ():
    """
    Convert pysam aligned segment Obj in pure python Obj, as pysam obj cannot be pickled
    """
    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

    def __init__(self, qname, qlen, qstart, qend, rname, rlen, rstart, rend, strand, align_len, mapq, align_score=None, **kwargs):
        """
        """
        self.qname = qname
        self.qlen = qlen
        self.qstart = qstart
        self.qend = qend
        self.rname = rname
        self.rlen = rlen
        self.rstart = rstart
        self.rend = rend
        self.strand = strand
        self.align_len = align_len
        self.mapq = mapq
        self.align_score = align_score

    def __repr__(self):
        return "Query:{}-{}:{} ({} pb) / Reference:{}-{}:{}({}) ({} pb) / Alignment len:{} / Mapq:{} / Align Score:{}".format(
            self.qname, self.qstart, self.qend, self.qlen,
            self.rname, self.rstart, self.rend, self.strand, self.rlen,
            self.align_len, self.mapq, self.align_score)
