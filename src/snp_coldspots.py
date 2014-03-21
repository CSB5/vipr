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
# smaller sizes might not give significant results depending on number
# of samples 50 works
#MIN_REGION_SIZE = 50
#MIN_REGION_SIZE = 2
#MIN_REGION_SIZE = 40

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
    parser.add_option("", "--minsize",
                      dest="min_size",
                      type="int",
                      help="Minimum region size (e.g. 40)")
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



def find_coldspot_regions(pooled_snps, seq_len, min_reg_size, excl_pos):
    """FIXME:add-doc-str
    """

    # otherwise math doesn't make sense
    snps_in_excl = []
    for s in pooled_snps:
        assert s.pos < seq_len
        if s.pos in excl_pos:
            snps_in_excl.append(s)
    if len(snps_in_excl):
        LOG.warn("%d of %d snps were in excl_pos (pos: %s)" % (
            len(snps_in_excl), len(pooled_snps), 
            ','.join([str(s.pos+1) for s in snps_in_excl])))
 
    all_snp_pos = set([s.pos for s in pooled_snps if s.pos not in excl_pos])
    
    # don't count twice (my conservative version)
    # snp_prob = len(all_snp_pos)/float(seq_len)
    # or
    # count all (as told by NN): len([s.pos for s in pooled_snps if s.pos not in excl_pos])
    snp_prob = len(all_snp_pos)/(float(seq_len)-len(excl_pos))
    # FIXME This ignores the fact that we can have 3 SNPs per pos.
    # Such cases will give probs>1 and the binom_test will fail
    #LOG.warn("snp_prob = %f" % snp_prob)

    coldspot_regions = []
    cur_region_start = 0

    # keep a list of positions where we stop and look back at the
    # stretch of SNP void regions. Those stops include excl_pos, which
    # only works because they are a range and can therefore never be a
    # coldspot on their own.
    stop_pos = [p for p in excl_pos if p<seq_len]
    stop_pos.extend(all_snp_pos)
    stop_pos = set(sorted(stop_pos))
    for cur_snp_pos in sorted(stop_pos):
        diff = cur_snp_pos - cur_region_start
        if diff >= min_reg_size:
            coldspot_regions.append((cur_region_start, cur_snp_pos-1))
            #print "DEBUG: %d - %d" % (cur_region_start, cur_snp_pos-1)
        cur_region_start = cur_snp_pos+1

    # needs trimming of min_reg_size above as well
    #diffs = []
    #for region in coldspot_regions:
    #    (region_start, region_end) = region
    #    diffs.append(region_end-region_start)
    #import numpy
    #print "diffs", numpy.mean(diffs), "3*", numpy.std(diffs)
    #import pdb; pdb.set_trace()

    bonf_fac = len(coldspot_regions)
    for region in coldspot_regions:
        (region_start, region_end) = region
        size = region_end-region_start+1
        # one tailed using reverse logic (does that make sense?)
        #pvalue = binom.sf(size-1, size, 1-snp_prob)
        # two tailed:
        try:
            # two sided pvalue_old = binom_test(0, size, snp_prob)
            # one sided:
            pvalue = binom.cdf(0, size, snp_prob)
        except:
            #import pdb; pdb.set_trace()
            raise
        #LOG.warn("coldspot (0 SNPs) region %d-%d (length %d): un-corrected pvalue = %g (bonf. factor is %g)" % (
        #        region_start, region_end, size, pvalue, bonf_fac))
        if pvalue * bonf_fac <= SIG_LEVEL:
            print "coldspot region %d-%d (length %d): bonferroni corrected pvalue = %g (uncorrected = %g)" % (
                region_start, region_end, size, pvalue * bonf_fac, pvalue)
        else:
            LOG.info("non-significant coldspot region %d-%d (length %d):"
                     " bonferroni corrected pvalue = %g (uncorrected = %g)" % (
                         region_start, region_end, size, pvalue * bonf_fac, pvalue))


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
            (opts.min_size, "minimum region size"),
            ]:
        if not required_opt:
            LOG.fatal("Missing %s argument" % opt_descr)
            sys.exit(1)
 
    LOG.info("Init sig-level=%f; min region size = %d" % (
        SIG_LEVEL, opts.min_size))

    excl_pos = []
    if opts.fexclude:
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Excluding %d positions" % len(excl_pos))

    snps = []
    for snp_file in snp_files:
        more_snps = snp.parse_snp_file(snp_file)
        snps.extend(more_snps)
        LOG.info("Parsed %d SNPs from %s" % (len(more_snps), snp_file))
        
    find_coldspot_regions(snps, opts.seq_len, opts.min_size, excl_pos)
    
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
