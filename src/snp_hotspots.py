#!/usr/bin/env python
"""Determine mutational hotspots
"""

#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser


# --- third-party imports
#
#from scipy.stats.distributions import binom
from scipy.stats import binom_test

#--- project specific imports
#
from lofreq import snp


SIG_LEVEL = 0.05


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


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    parser.add_option("", "--winsize",
                      dest="win_size",
                      type="int",
                      help="Window size to use (e.g. 20)")
    parser.add_option("", "--step",
                      dest="step_size",
                      type="int",
                      help="Step size to use (e.g. 5)")
    parser.add_option("", "--snpfile",
                      dest="snp_file",
                      help="SNP file to read")
    parser.add_option("", "--seqlen",
                      dest="seq_len",
                      type="int",
                      help="Sequence length")
    return parser



def find_hotspot_windows(snps, seq_len, win_size, step_size):
    """FIXME
    """

    
    bonf_fac = seq_len/float(win_size)
    snp_prob = len(snps)/float(seq_len)

    LOG.info("Using a window size of %d with a step size of %d. SNP prob is %g" % (
            win_size, step_size, snp_prob))
    
    curpos = 0
    #for win_start in range(0, seq_len, win_size):
    while curpos < seq_len-win_size:
        win_start = curpos
        num_snps_in_win = len([s for s in snps
                               if s.pos>=win_start and s.pos<=(win_start+win_size)])
        if num_snps_in_win > 0:
            #print win_start, win_start+win_size, num_snps_in_win

            if num_snps_in_win > win_size:
                LOG.warn("Hack to make num_snps_in_win>win_size possible")
                num_snps_in_win = win_size
            # one tailed: reports pos with even 1 at <0.05 on den2 replicates
            # pvalue = binom.sf(1, win_size, len(snps)/float(SEQ_LEN))
            # two tailed:
            try:
                pvalue = binom_test(num_snps_in_win, win_size, snp_prob)
            except ValueError:
                LOG.fatal("The following failed: binom_test(num_snps_in_win=%d, win_size=%d, snp_prob=%d)" % (
                    num_snps_in_win, win_size, snp_prob))
                raise
            #LOG.warn("DEBUG: window %d-%d has %d SNPs which translates into a raw pvalue of %g" % (
            #        win_start, win_start+win_size, num_snps_in_win, pvalue))
            if pvalue * bonf_fac <= SIG_LEVEL:
                print "window %d-%d has %d SNPs which translates into a bonferroni corrected pvalue of %g" % (
                    win_start, win_start+win_size, num_snps_in_win, pvalue * bonf_fac)
                
        curpos += step_size


def main():
    """
    The main function
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

    for (required_opt, opt_descr) in [(opts.win_size, "window size"), 
                                      (opts.step_size, "step size"),
                                      (opts.snp_file, "SNP file"),
                                      (opts.seq_len, "sequence length")]:
        if not required_opt:
            LOG.fatal("Missing %s argument" % opt_descr)
            sys.exit(1)

    snps = snp.parse_snp_file(opts.snp_file)
    LOG.info("Parsed %d SNPs from %s" % (len(snps), opts.snp_file))

    find_hotspot_windows(snps, opts.seq_len, opts.win_size, opts.step_size)
    
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
