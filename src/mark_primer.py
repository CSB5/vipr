#!/usr/bin/env python
"""The script will take a BAM file and fw/rv primer positions and
remove any reads perfectly overlapping these primer positions. It will
also clip bases partially overlapping with the primer regions.
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, OptionGroup
from collections import namedtuple


#--- third-party imports
#
import pysam
import cigar
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


DEFAULT_PRIMER_LEN = 25


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
    parser.add_option("", "--force",
                      action="store_true",
                      dest="force_overwrite",
                      help="Optional: force overwriting of files")
    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="Optional: enable debugging")
    parser.add_option("-i", "--bam-src",
                      dest="fbam_in",
                      help="BAM input file")
    parser.add_option("-o", "--bam-dst",
                      dest="fbam_out",
                      help="BAM output file")
    parser.add_option("-p", "--primer",
                      dest="primer_pos_file",
                      help="Primer start position file Format:"
                      " 1-based-start-pos orientation-F-or-R.")
    parser.add_option("", "--primer-len",
                      dest="primer_len",
                      type="int",
                      default=DEFAULT_PRIMER_LEN,
                      help="Don't ignore all of peak positions"
                      " but only the actual primer pos, i.e. its"
                      " corresponding beginning until primer length"
                      " (default %d)" % (DEFAULT_PRIMER_LEN))

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



def print_alignedread_info(alignedread):
    """Debugging only"""

    print "- pos, aend-1, alen = %d %d %s" % (
        alignedread.pos, alignedread.aend-1, alignedread.alen)
    print "- cigar = %s" % (alignedread.cigar)
    print "- decoded cigar = %s" % (
        ''.join(list(cigar.decoded_ops(alignedread.cigar))))
    try:
        alnlen = cigar.aligned_length(alignedread.cigar)
    except ValueError:
        alnlen =  -1
        print "- cigar aln-len = %d" % (alnlen)
        print "- read sequence = %s" % (alignedread.query)
        print "- aligned parts = %s" % (alignedread.query)


        
def clip_primer_overlap_bases(alignedread, primer_overlap_peak_pos, primer_len):
    """Clip alignread's bases overlapping with primer pos by changing
    cigar string
    """

    if alignedread.is_reverse:
        # WARNING: partial code duplication with else branch
        #
        # this is a reverse read overlapping with it's 5p end,
        # so starting from there clip bases, skipping any
        # non-match ops. do this p+primer_len -
        # alignedread.pos + 1 times
        num_clip = primer_overlap_peak_pos + primer_len - alignedread.pos + 1
        LOG.debug("rv read starting at %d overlapping with fw peak start %d+%d by %d: %s" % (
            alignedread.pos, primer_overlap_peak_pos, primer_len, num_clip, alignedread))
        new_cigar_decoded = []
        for op in cigar.decoded_ops(alignedread.cigar):
            if num_clip == 0 or op == 'D':
                new_cigar_decoded.append(op)
            else:
                # replace everything else with S
                new_cigar_decoded.append('S')
                if op in 'MI=X':
                    num_clip -= 1

    else:
        # WARNING: partial code duplication with if branch
        #
        # this is a forward read overlapping with its' 3p end, so
        # starting from there clip bases, skipping any non-match
        # ops. do this p+primer_len - alignedread.pos + 1 times
        num_clip = alignedread.aend - primer_overlap_peak_pos + primer_len + 1
        LOG.debug("fw read ending at %d overlapping with rv peak end %d-%d by %d: %s" % (
            alignedread.aend, primer_overlap_peak_pos, primer_len, num_clip, alignedread))
        new_cigar_decoded = []
        for op in list(cigar.decoded_ops(alignedread.cigar))[::-1]:
            if num_clip == 0 or op == 'D':
                new_cigar_decoded.insert(0, op)
            else:
                # replace everything else with S
                new_cigar_decoded.insert(0, 'S')
                if op in 'MI=X':
                    num_clip -= 1
        
    new_cigar = cigar.parse(
        cigar.cigar_from_decoded_ops(new_cigar_decoded))
    # According to spec http://samtools.sourceforge.net/SAM1.pdf:
    # "Sum of lengths of the M/I/S/=/X operations shall equal the
    # length of SEQ"
    cigar_len = sum([1 for op in new_cigar_decoded if op in 'MISX='])
    assert cigar_len == alignedread.rlen, (
        "read length derived from new cigar (%s -> %d) mismatches rlen=%d for %s read: %s" % (
            new_cigar, cigar_len, alignedread.rlen,
            "reverse" if alignedread.is_reverse else "forward", alignedread))

    alignedread.cigar = new_cigar

    

def mark_primer(sam_in, sam_out,
                peaks_fw_start_pos, peaks_rv_end_pos, primer_len,
                plpregion=None):
    """FIXME:add-doc
    """

    plpref, plpstart, plpend = plpregion

    num_reads_out = 0
    num_reads_in = 0
    num_fw_peak_start = 0
    num_rv_peak_end = 0 # i.e. right-most rv peak pos
    num_fw_primer_overlap = 0
    num_rv_primer_overlap = 0

    # read arrival rates using left-most coordiante of fw reads and
    # the right-most coordinate for rv reads
    fw_read_start_rate = dict()
    rv_read_end_rate = dict()

    ctr = 0
    for alignedread in sam_in.fetch(plpref, plpstart, plpend):
        num_reads_in += 1
        ctr += 1
        if ctr%500000 == 0:
            LOG.info("Processing read no %d..." % ctr)
            
        # Most important data structure:
        # http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.AlignedRead
        #
        # refseq[alignedread.pos:alignedreas.aend] will give you the
        # aligned part of the refseq. this also means aend should not
        # be used as position as it's the open end of a slice (i.e.
        # use -1 where necessary)


        # debugging
        #
        if False:
            if alignedread.cigar[0][0] != 0 or alignedread.cigar[-1][0] != 0:
                print "Leading/trailing non-match in read %s" % (alignedread)
            elif len(alignedread.cigar)>1:# or alignedread.cigar[0][0] != 0:
                print "Imperfect mapping of %s" % (alignedread)
            print_alignedread_info(alignedread)
            #import pdb; pdb.set_trace()
 

        # if unmapped, write and skip to next read
        #
        if alignedread.is_unmapped:
            LOG.debug("Leaving unmapped read untouched: %s" % alignedread)
            sam_out.write(alignedread)
            num_reads_out += 1
            continue


        # easy cases first: mark as duplicate if
        # 1. reverse and read rightmost coordinate is
        # in rv-primer-end or
        # 2. forward and leftmost coordinate is
        # in fw-primer-start
        #
        # if any of this is true, mark as pcr duplicate, output and
        # continue with next read
        #
        if alignedread.is_reverse:
            rv_read_end_rate[alignedread.aend-1] = rv_read_end_rate.get(
                alignedread.aend-1, 0) + 1
            
            if alignedread.aend-1 in peaks_rv_end_pos:
                num_rv_peak_end += 1
                LOG.debug("Marking reverse read ending at rv peak end"
                          " as duplicate: %s" % (alignedread))
                alignedread.flag |= 0x400
                sam_out.write(alignedread)
                num_reads_out += 1
                continue
        else:
            fw_read_start_rate[alignedread.pos] = fw_read_start_rate.get(
                alignedread.pos, 0) + 1
            
            if alignedread.pos in peaks_fw_start_pos:
                num_fw_peak_start += 1
                LOG.debug("Marking forward read starting at fw peak start"
                          " as duplicate: %s" % (alignedread))
                alignedread.flag |= 0x400
                sam_out.write(alignedread)
                num_reads_out += 1
                continue

        
        # the more complicated cases: if read is on the opposite site of a peak
        # and reads into the primer region, then clip its bases
        #
        #
        if alignedread.is_reverse: # = alignedread.flag & 0x10:
            primer_overlap_peak_pos = [p for p in peaks_fw_start_pos
                                       if alignedread.pos >= p
                                       and alignedread.pos <= p+primer_len]
            num_fw_primer_overlap += len(primer_overlap_peak_pos)
        else:
            primer_overlap_peak_pos = [p for p in peaks_rv_end_pos
                                       if alignedread.aend-1 <= p
                                       and alignedread.aend-1 >= p-primer_len]
            num_rv_primer_overlap += len(primer_overlap_peak_pos)

        # if no overlap is found, then it's a normal read: output and
        # continue
        #
        if not len(primer_overlap_peak_pos):
            if False:
                LOG.debug("Leaving normal read untouched: %s" % (alignedread))
            sam_out.write(alignedread)
            num_reads_out += 1
            continue
        
        # LOG.critical("Skipping clipping"); continue
        
        for pos in primer_overlap_peak_pos:
            #oldcigar = alignedread.cigar
            clip_primer_overlap_bases(alignedread, pos, primer_len)
            #if alignedread.cigar == oldcigar:
            #    LOG.critical("No change in cigar after clipping. Might be ok if everything's clipped already")
            #    import pdb; pdb.set_trace()                
        sam_out.write(alignedread)
        num_reads_out += 1

    if num_reads_in == 0:
        LOG.fatal("No reads in input!")
        sys.exit(1)
    if num_reads_out == 0:
        LOG.critical("No reads passed!")
        
    LOG.info("# Fw-reads starting at fw-peak start, therefore marked as dup = %d" % num_fw_peak_start)
    LOG.info("# Rv-reads ending at rv-peak end, therefore marked as dup = %d" % num_rv_peak_end)
    LOG.info("# Clipped rv-reads with overlap in fw-primer = %d" % num_fw_primer_overlap)
    LOG.info("# Clipped fw-reads with overlap in rv-primer = %d" % num_rv_primer_overlap)


    for (rate_dict, name) in [(fw_read_start_rate, "fw start rate"),
                              (rv_read_end_rate, "rv end rate")]:
        mean_val = mean(rate_dict.values())
        std_val = std(rate_dict.values())
        thresh = mean_val + 3*std_val
        # or 5% quartile instead of top 3*sigma?
        outlier_pos = sorted([p for (p, v) in rate_dict.iteritems() if v>thresh])
        for pos in outlier_pos:
            LOG.info("Bonus: (re-)detected peak in %s at pos %d with rate %d" % (
                name, pos+1, rate_dict[pos]))


            
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
            (opts.fbam_in, 'BAM input file', 'in', True),
            (opts.fbam_out, 'BAM output file', 'out', True),
            (opts.primer_pos_file, 'Primer positions file', 'in', True),
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


    # primer length
    #
    if opts.primer_len < 0:
        LOG.fatal("Negative primer length does not make sense")
        sys.exit(1)
    primer_len = opts.primer_len


    # parse peaks and prep
    #
    peaks = pp.parse_primer_pos(open(opts.primer_pos_file))
    # summarize fw-start and rv-end positions in list. actually don't
    # need anything else (i.e. the end positions)
    peaks_fw_start_pos = [p.pos for p in peaks if p.ori == 'F']
    peaks_rv_end_pos = [p.pos for p in peaks if p.ori == 'R']

    sam_in = pysam.Samfile(opts.fbam_in, "rb")
    sam_out = pysam.Samfile(opts.fbam_out, "wb", template=sam_in)

    mark_primer(sam_in, sam_out,
                peaks_fw_start_pos, peaks_rv_end_pos, primer_len,
                (plpref, plpstart, plpend))

    if sam_in != sys.stdin:
        sam_in.close()
    if sam_out != sys.stdout:
        sam_out.close()



if __name__ == "__main__":

    pysam_version = tuple(int(x) for x in pysam.__version__.split('.'))
    if pysam_version < (0 , 6):
        LOG.warning("Using untested pysam version (%s)." % (pysam_version))

    main()
    LOG.info("Successful exit")



