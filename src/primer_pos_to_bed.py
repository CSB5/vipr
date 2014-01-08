#!/usr/bin/env python
"""This script will mask primer regions in a sequence.
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
from itertools import count, groupby


#--- third-party imports
#
from lofreq import sam

#--- project specific imports
#
import primer_pos as pp

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013 Genome Institute of Singapore"
__license__ = "GPL2"
__credits__ = [""]
__status__ = "eternal alpha"


DEFAULT_PRIMER_LEN = 25


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def lists_to_ranges(lists):
    """Inspired by http://stackoverflow.com/questions/3429510/pythonic-way-to-convert-a-list-of-integers-into-a-string-of-comma-separated-rang/3430231#3430231

    Holy cow...
    """

    return (list(x) for _, x in groupby(lists, lambda x, c=count(): next(c)-x))
     
    
def primer_positions_to_incl_bed(primer_positions, bed_fh, primer_len, 
                                 seq_len, seq_name):
    """FIXME:add-doc"""
    assert primer_len > 0 and seq_len > 0

    # create a dict that with keys ranging from 0 to seq_len-1. value
    # will be one if positions is within primer pos. default to 0
    match_dict = dict(zip(range(seq_len), seq_len*[0]))
    
    for primer_pos in primer_positions:
        if primer_pos.ori == 'F':
            for p in range(primer_pos.pos, primer_pos.pos+primer_len):
                match_dict[p] = 1
        elif primer_pos.ori == 'R':
            for p in range(primer_pos.pos, primer_pos.pos-primer_len, -1):
                match_dict[p] = 1           
        else:
            raise ValueError

    for r in lists_to_ranges([p for p in sorted(match_dict.iterkeys())
                              if match_dict[p]==0]):
        bed_fh.write("%s\t%d\t%d\n" % (seq_name, r[0], r[-1]+1))

        
def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog:\n" \
        + __doc__ + "\n" \
        "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="Optional: be verbose")
    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="Optional: enable debugging")
    parser.add_option("", "--force",
                      dest="force_overwrite",
                      action="store_true",
                      help="Force overwriting of files")
    parser.add_option("-i", "--primer-pos",
                      dest="primer_pos_file",
                      help="Primer position file")
    parser.add_option("-p", "--primer-len",
                      dest="primer_len",
                      default=DEFAULT_PRIMER_LEN,
                      type="int",
                      help="Primer length (default %d)" % DEFAULT_PRIMER_LEN)
    parser.add_option("-b", "--bam",
                      dest="bam_file",
                      help="Corresponding BAM file to derive --seqname and --seqlen automatically")
    parser.add_option("-o", "--bed",
                      dest="bed_file",
                      default="-",
                      help="Convert to this bed-file listing include regions"
                      " (use for e.g. samtools mpileup -l). Use '-' for stdout (default)")
    parser.add_option("", "--seqlen",
                      type="int",
                      dest="seq_len",
                      help="Sequence length")
    parser.add_option("", "--seqname",
                      dest="seq_name",
                      help="Sequence name")

    return parser


            
def main():
    """The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()
    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    # file check
    #
    for (filename, descr, direction, mandatory) in [
            (opts.bam_file, 'BAM input file', 'in', False),
            (opts.primer_pos_file, 'primer position input file', 'in', True),
            (opts.bed_file, 'bed output file', 'out', False),
            ]:

        if not mandatory and not filename:
            continue

        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)

        if filename == '-':
            continue

        if direction == 'in' and not os.path.exists(filename):
            LOG.fatal(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)

        if direction == 'out' and os.path.exists(filename) and not opts.force_overwrite:
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)

    # primer length
    #
    if opts.primer_len < 0:
        LOG.fatal("Negative primer length does not make sense")
        sys.exit(1)

    if not opts.bam_file:
        if not opts.seq_len or not opts.seq_name:
            LOG.fatal("Missing sequence name or sequence length argument")
            sys.exit(1)
        seq_len = opts.seq_len
        seq_name = opts.seq_name
    else:
        if opts.seq_len or opts.seq_name:
            LOG.fatal("BAM file given, so will derive seqlen and"
                      " seqname automatically, which were however"
                      " also given as arguments")
            sys.exit(1)
        sam_header = sam.sam_header(opts.bam_file)
        sq_list = sam.sq_list_from_header(sam_header)
        assert len(sq_list)==1, (
            "Can only work with one sequence but found %d in %s" % (
                (len(sq_list), opts.bam_file)))
        seq_name =  sq_list[0]
        seq_len = sam.len_for_sq(sam_header, seq_name)
        
    if opts.bed_file == "-":
        bed_fh = sys.stdout
    else:
        bed_fh = open(opts.bed_file, 'w')

    
    primer_positions = pp.parse_primer_pos(open(opts.primer_pos_file))
    primer_positions_to_incl_bed(
        primer_positions, bed_fh, opts.primer_len,
        seq_len, seq_name)
    
    if bed_fh != sys.stdout:
        bed_fh.close()
        
if __name__ == "__main__":
    main()
    LOG.info("Successful exit")



