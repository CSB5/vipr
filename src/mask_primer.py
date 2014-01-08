#!/usr/bin/env python
"""This script will mask primer regions in a sequence
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser



#--- third-party imports
#
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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


def mask_primer_in_seq(seqrec,
                       primers_fw_start_pos, primers_rv_end_pos, 
                       primer_len, ign_index_error=False):
    """Will mask primer regions by modifing seqrec such that each
    region given by fw_peak+primer_len-1 and rv_peak-primer_len-1 is
    turned into Ns. If regions are outside of the sequence an
    IndexError will be raised unless ign_index_error=True

    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> seqrec = SeqRecord(Seq("ACA", generic_dna))
    >>> mask_primer_in_seq(seqrec, [1], [], 1)
    >>> print seqrec.seq
    ANA
    >>> seqrec = SeqRecord(Seq("ACA", generic_dna))
    >>> mask_primer_in_seq(seqrec, [], [1], 1)
    >>> print seqrec.seq
    ANA
    >>> seqrec = SeqRecord(Seq("ACGTACGTACGT", generic_dna))
    >>> mask_primer_in_seq(seqrec, [2], [9], 2)
    >>> print seqrec.seq
    ACNNACGTNNGT
    >>> mask_primer_in_seq(seqrec, [], [100], 1, True)
    >>> print seqrec.seq
    ACNNACGTNNGT
    >>> mask_primer_in_seq(seqrec, [], [100], 1)
    Traceback (most recent call last):
      ...
    IndexError: list assignment index out of range
    
    See http://stackoverflow.com/questions/12592/can-you-check-that-an-exception-is-thrown-with-doctest-in-python
    """

    seq = list(seqrec.seq)
    for p in primers_fw_start_pos:
        assert primers_fw_start_pos >= 0
        for i in range(primer_len):
            try:
                seq[p+i] = 'N'
            except IndexError:
                # okay if we are running over the edge but start was ok
                if p >= len(seq):
                    raise
                else:
                    break
    for p in primers_rv_end_pos:
        for i in range(primer_len):
            try:
                seq[p-i] = 'N'
            except IndexError:
                # okay if we are running over the edge but start was ok
                if p < 0:
                    raise
                else:
                    break
    seqrec.seq = Seq(''.join(seq), generic_dna)


def mask_primer(fh_fa_in, fh_fa_out,
                primers_fw_start_pos, primers_rv_end_pos, primer_len,
                in_fmt="fasta", out_fmt="fasta"):
    """Call mask_primer_in_seq for each sequence record in file hande
    fh_fa_in and report masked seq to fh_fa_out
    """

    for seqrec in SeqIO.parse(fh_fa_in, in_fmt):
        mask_primer_in_seq(seqrec,
                           primers_fw_start_pos, primers_rv_end_pos,
                           primer_len)
        SeqIO.write([seqrec], fh_fa_out, out_fmt)

        
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
    parser.add_option("-i", "--fa-src",
                      dest="f_fa_in",
                      help="Input fasta file")
    parser.add_option("-o", "--fa-dst",
                      dest="f_fa_out",
                      help="Output fasta file")
    parser.add_option("-p", "--primers",
                      dest="primers",
                      help="Primer positions file, i.e. coverage peaks"
                      " due to primers. Format:"
                      " 0-based-start-pos include-end-pos direction-F-or-R.")
    parser.add_option("", "--primer-len",
                      dest="primer_len",
                      type="int",
                      default=DEFAULT_PRIMER_LEN,
                      help="Don't ignore all of peak positions"
                      " but only the actual primer pos, i.e. its"
                      " corresponding beginning until primer length"
                      " (default %d)" % (DEFAULT_PRIMER_LEN))
    parser.add_option("", "--force",
                      dest="force_overwrite",
                      action="store_true",
                      help="Force overwriting of files")
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
            (opts.f_fa_in, 'Fasta input file', 'in', True),
            (opts.f_fa_out, 'Fasta output file', 'out', True),
            (opts.primers, 'Primer position file', 'in', True),
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

        if direction == 'out' and os.path.exists(filename) \
          and not opts.force_overwrite:
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)

    # primer length
    #
    if opts.primer_len < 0:
        LOG.fatal("Negative primer length does not make sense")
        sys.exit(1)
    primer_len = opts.primer_len


    # parse primers and prep
    #
    primer_pos = pp.parse_primer_pos(open(opts.primers))
    # summarize fw-start and rv-end positions in list. actually don't
    # need anything else (i.e. the end positions)
    primers_fw_start_pos = [p.pos for p in primer_pos if p.ori == 'F']
    primers_rv_end_pos = [p.pos for p in primer_pos if p.ori == 'R']
                
    if opts.f_fa_in == '-':
        fh_fa_in = sys.stdin
    else:
        fh_fa_in = open(opts.f_fa_in, 'r')
    if opts.f_fa_out == '-':
        fh_fa_out = sys.stdout
    else:
        fh_fa_out =  open(opts.f_fa_out, 'w')

    mask_primer(fh_fa_in, fh_fa_out,
                primers_fw_start_pos, primers_rv_end_pos, primer_len)

    if fh_fa_in != sys.stdin:
        fh_fa_in.close()
    if fh_fa_out != sys.stdout:
        fh_fa_out.close()



if __name__ == "__main__":
    main()
    LOG.info("Successful exit")



