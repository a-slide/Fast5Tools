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

    def add_read (self, pysam_read):
        """"""
        self.read_list.append (Read (pysam_read))
        self.nread += 1

    def best_read (self, based_on = "score"):
        """"""
        if self.nread == 0:
            return None

        if self.nread == 1:
            return self.read_list [0]

        if based_on == "score":
            return sorted(self.read_list, key=lambda x: x.score, reverse=True)

        elif based_on == "mapq":
            return sorted(self.read_list, key=lambda x: x.mapq, reverse=True)

        elif based_on == "rlen":
            return sorted(self.read_list, key=lambda x: x.rlen, reverse=True)

    def __repr__(self):
        """ Readable description of the object """
        m = ""
        for r in self.read_list:
            m +="\t\t{}\n".format(r)
        return (m)

class Read ():
    """
    Convert pysam aligned segment Obj in pure python Obj, as pysam obj cannot be pickled
    """
    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

    def __init__(self, pysam_read, **kwargs):
        """ Store initial values and init counters
        """
        try:
            self.qname = pysam_read.query_name
            self.rname = pysam_read.reference_name
            self.start = int (pysam_read.reference_start)+1
            self.end = int (pysam_read.reference_end)
            self.strand = "-" if pysam_read.is_reverse else "+"
            self.rlen = int (pysam_read.infer_read_length())
            self.mapq = int (pysam_read.mapping_quality)
            self.flag = pysam_read.flag
            self.score = int(pysam_read.get_tag("AS"))

        except TypeError as E:
            print (pysam_read)
            print (E)

    def __repr__(self):
        return "Reference:{}-{}:{}({})  Read length:{}  Mapq:{},  Score:{}  Flag:{}".format(
            self.rname, self.start, self.end, self.strand, self.rlen, self.mapq, self.score, self.flag)
