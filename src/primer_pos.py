#!/usr/bin/env python
"""FIXM:add-doc"""

#--- standard library imports
#
from collections import namedtuple
import logging


#--- third-party imports
#
#/

#--- project specific imports
#
#/

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013 Genome Institute of Singapore"
__license__ = "GPL2"
__credits__ = [""]
__status__ = "eternal alpha"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


PrimerPos = namedtuple('PrimerPos', ['pos', 'ori'])


def write_primer_pos(primer_pos, fh):
    """FIXME:add-doc"""

    # header
    fh.write("# Primer positions and orientation are reported per line as:\n")
    fh.write("# <pos> <ori>\n")
    fh.write("# where <pos> is the start-position (inclusive, one-based) and\n")
    fh.write("# <ori> is one of [FR], i.e. forward or reverse\n")
    fh.write("# Note, that F primers extend to <pos>+ and R primers to <pos>-\n")
    
    for p in primer_pos:
        fh.write("%d %c\n" %(p.pos+1, p.ori))# internally: zero-offset
        
        
def parse_primer_pos(fh):
    """FIXME:add-doc"""
    
    primer_pos = []    
    for line in fh:
        line = line.rstrip()
        if len(line)==0 or line.startswith("#"):
            continue

        line_split = line.split()
        assert len(line_split)==2, (
            "Expected exactly two values but got: '%s'" % (line_split))
        (pos, ori) =  line_split
        pos = int(pos)-1# internally: zero-offset
        assert pos >= 0 and ori in 'FR', (
            "Oups...didn't understand the following line: %s" % line)
        primer_pos.append(PrimerPos(pos, ori))
    return primer_pos

        


if __name__ == "__main__": 
    import sys
    sys.stderr.write("Go away. Nothing to see here")

