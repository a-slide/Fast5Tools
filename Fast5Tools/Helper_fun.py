# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def stdout_print (*args):
    """reproduce print with stdout.write
    """
    sys.stdout.write(" ".join(str(a) for a in args))
    sys.stdout.flush()
