#!/usr/bin/env python
"""Determine mutational coldspots from pooled SNP files
"""

#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser


# --- third-party imports
#
# /
#from scipy.stats.distributions import binom
from scipy.stats import binom_test
from scipy.stats import binom

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
    parser.add_option("", "--seqlen",
                      dest="seq_len",
                      type="int",
                      help="Sequence length")
    parser.add_option("", "--exclude",
                      dest="fexclude",
                      help="Optional: exclude positions files"
                      " (zero-based, half-open; like bed without chromosome)")
    return parser




def read_exclude_pos_file(fexclude):
    """Taken from lofreq_snpcaller.py

    Parse file containing ranges of positions to exclude and return
    positions as list. File format is like bed without chromosome
    """

    excl_pos = []
    fhandle = open(fexclude, 'r')
    for line in fhandle:
        if line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue

        start = int(line.split()[0])
        end = int(line.split()[1])
        # ignore the rest

        assert start < end, ("Invalid position found in %s" % fexclude)

        excl_pos.extend(range(start, end))
    fhandle.close()

    return excl_pos



def find_coldspot_regions(pooled_snps, seq_len, excl_pos):
    """FIXME:add-doc-str
    """
        
    # remove snvs in excl pos
    #
    snps_in_excl = []
    for s in pooled_snps:
        assert s.pos < seq_len
        if s.pos in excl_pos:
            snps_in_excl.append(s)
    if len(snps_in_excl):
        LOG.warn("%d of %d snps were in excl_pos (pos: %s)" % (
            len(snps_in_excl), len(pooled_snps), 
            ','.join([str(s.pos+1) for s in snps_in_excl])))
 
    
    # WARN question is whether to count SNVs occuring at the same pos
    # as one len(set(non_excl_pos)) or multiple times (no set()).
    # NN says to count all
    #
    # WARN This ignores the fact that we can have 3 SNPs per pos.
    # Such cases might give probs>1 and the binom_test will fail
    #LOG.warn("snp_prob = %f" % snp_prob)
    #
    non_uniq_incl_snp_pos = [s.pos for s in pooled_snps if s.pos not in excl_pos]
    if False:
        snp_prob = len(set(non_uniq_incl_snp_pos))/(float(seq_len)-len(excl_pos))
    else:
        snp_prob = len(non_uniq_incl_snp_pos)/(float(seq_len)-len(excl_pos))
    incl_snp_pos = set(non_uniq_incl_snp_pos)

    
    # compute coldspot region candidates
    #
    # keep a list of positions where we stop and look back at the
    # stretch of SNP void regions. Those stops include excl_pos, which
    # only works because they are a range and can therefore never be a
    # coldspot on their own.
    coldspots_candidates = []
    cur_region_start = 0
    stop_pos = [p for p in excl_pos if p<seq_len]
    stop_pos.extend(incl_snp_pos)
    stop_pos = set(sorted(stop_pos))
    for cur_snp_pos in sorted(stop_pos):
        if cur_region_start+1<cur_snp_pos:
            coldspots_candidates.append((cur_region_start, cur_snp_pos-1))
        cur_region_start = cur_snp_pos+1

        
    bonf = len(coldspots_candidates)

    
    # compute minimum length of snv-void region considered significant
    #
    min_len = 1
    while 1:
        # 0 = zero successes/snvs observerd
        pv = binom.cdf(0, min_len, snp_prob)
        if pv*bonf < 0.05:
            break
        min_len += 1
    LOG.info("minimum required length = %d" % min_len)

    
    print "#length\tcorr. pv (%g)\tfrom\tto"
    for region in coldspots_candidates:
        (region_start, region_end) = region
        size = region_end-region_start+1
        try:
            # two sided pvalue_old = binom_test(0, size, snp_prob)
            # one sided:
            pvalue = binom.cdf(0, size, snp_prob)
        except:
            #import pdb; pdb.set_trace()
            raise
        if size<min_len:
            cat = "REJECT"
        else:
            cat = "PASSED"
        LOG.info("coldspot region %d-%d (length %d): newly computed and bonferroni corrected pvalue = %g" % (
                 region_start, region_end, size, pvalue * bonf))
        if cat == "PASSED":
            print "%d\t%g\t%d\t%d" % (size, pvalue * bonf, region_start, region_end)


def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    snp_files = args
    if len(snp_files)<2:
        parser.error("Need more than two SNP files.")
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    for (required_opt, opt_descr) in [
            (opts.seq_len, "sequence length"),
            ]:
        if not required_opt:
            LOG.fatal("Missing %s argument" % opt_descr)
            sys.exit(1)
 
    LOG.info("Init sig-level=%f" % (SIG_LEVEL))

    excl_pos = []
    if opts.fexclude:
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Excluding %d positions" % len(excl_pos))

    snps = []
    for snp_file in snp_files:
        more_snps = snp.parse_snp_file(snp_file)
        snps.extend(more_snps)
        LOG.info("Parsed %d SNPs from %s" % (len(more_snps), snp_file))
        
    print "#coldspots"
    print "#excluding positions listed in %s" % opts.fexclude
    print "#considering snvs listed in %s" % (', '.join(snp_files))
    find_coldspot_regions(snps, opts.seq_len, excl_pos)
    
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
