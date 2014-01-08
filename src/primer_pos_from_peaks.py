#!/usr/bin/env python
"""The script will determine primer positions (coverage peaks) on the
forward and reverse strand of a BAM file. This done by first
extracting the read start/end rate for each strand and then applying a
simple 3-sigma cutoff.

NOTES:
- This code was written for single, short (viral) genomes. Your milage
on multichromosome and larger genomes will vary.
- In freak cases with indels in your mapping at primer positions, the
actual length of the reported peak might not necessarily be equal to
your read-length
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, OptionGroup


#--- third-party imports
#
import pysam
from numpy import mean, std

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
    parser.add_option("-i", "--bam",
                      dest="bam_file",
                      help="BAM input file")
    parser.add_option("-o", "--pos",
                      dest="primer_pos_file",
                      help="Primer positions output file")

    region_group = OptionGroup(parser, "Region Options", "")
    region_group.add_option("", "--ref",
                      dest="plpref",
                      help="Optional:Pileup reference (chromosome/sequence)."
                      " Only need if pileup start and/or end are set")
    region_group.add_option("", "--start",
                      dest="plpstart",
                      type="int",
                      help="Optional:Pileup start (if set, needs ref as well)")
    region_group.add_option("", "--end",
                      dest="plpend",
                      type="int",
                      help="Optional:Pileup end  (if set, needs ref as well)")
    parser.add_option_group(region_group)

    return parser

    

def find_primer_positions(sam, plpregion=None):
    """FIXME:add-doc
    """

    plpref, plpstart, plpend = plpregion
    # read arrival rates using left-most coordiante of fw reads and
    # the right-most coordinate for rv reads
    fw_read_start_rate = dict()
    rv_read_end_rate = dict()

    ctr = 0
    for alignedread in sam.fetch(plpref, plpstart, plpend):
        ctr += 1
        if ctr % 500000 == 0:
            LOG.info("Processing read no %d..." % ctr)
            
        # Most important data structure:
        # http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.AlignedRead
        #
        # refseq[alignedread.pos:alignedreas.aend] will give you the
        # aligned part of the refseq. this also means aend should not
        # be used as position as it's the open end of a slice (i.e.
        # use -1 where necessary)

        if alignedread.is_unmapped:
            continue
        
        if alignedread.is_reverse:
            # aend is the end of a slice. we want the pos. therefore -1
            rv_read_end_rate[alignedread.aend-1] = rv_read_end_rate.get(
                alignedread.aend-1, 0) + 1
        else:
            fw_read_start_rate[alignedread.pos] = fw_read_start_rate.get(
                alignedread.pos, 0) + 1

    primer_positions = []
    for (rate_dict, ori) in [(fw_read_start_rate, "F"),
                              (rv_read_end_rate, "R")]:
        mean_val = mean(rate_dict.values())
        std_val = std(rate_dict.values())
        thresh = mean_val + 3*std_val
        LOG.info("Strand %c threshold: %f = %f * 3*%f" % (
            ori, thresh, mean_val, std_val))
        for (pos, val) in rate_dict.iteritems():
             if val<=thresh:
                 continue
             primer_positions.append(pp.PrimerPos(pos, ori))
             LOG.info("%c pos %d with rate %d exceeds threshold" % (
                 ori, pos, val))
    return sorted(primer_positions, key = lambda x: x.pos)


            
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
            (opts.bam_file, 'BAM input file', 'in', True),
            (opts.primer_pos_file, 'Primer positions output file', 'out', True),
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

        if direction == 'out' and os.path.exists(filename):
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)

    # region arguments
    #
    plpref = None
    plpstart = None
    plpend = None
    if opts.plpstart or opts.plpend or opts.plpref:
        if None in [opts.plpstart, opts.plpend, opts.plpref]:
            LOG.fatal("If one of pileup-ref, -start or -end is given,"
                      " then all three must be given")
            sys.exit(1)
        if opts.plpstart < 1 or opts.plpstart >= opts.plpend:
            LOG.fatal("Pileup start and end coordinates don't make sense")
            sys.exit(1)

        plpref = opts.plpref
        plpstart = opts.plpstart-1
        plpend = opts.plpend-1

        
    bam_fh = pysam.Samfile(opts.bam_file, "rb") 
    if opts.primer_pos_file == '-':
        primer_pos_fh = sys.stdout
    else:
        primer_pos_fh = open(opts.primer_pos_file, 'w')

    primer_positions = find_primer_positions(
        bam_fh, (plpref, plpstart, plpend))

    pp.write_primer_pos(primer_positions, primer_pos_fh)
    
    if primer_pos_fh != sys.stdout:
        primer_pos_fh.close()
    if bam_fh != sys.stdin:
        bam_fh.close()



if __name__ == "__main__":

    pysam_version = tuple(int(x) for x in pysam.__version__.split('.'))
    if pysam_version < (0 , 6):
        LOG.warning("Using untested pysam version (%s)." % (pysam_version))

    main()
    LOG.info("Successful exit")



