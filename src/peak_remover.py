#!/usr/bin/env python
"""Removes reads matching a certain start position from BAM file. You
can in theory achieve the same with a bamtools json filter.

See also peak_finder.py
"""

#--- standard library imports
#
from __future__ import division
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import os
import subprocess
from collections import namedtuple

#--- third-party imports
#

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013 Genome Institute of Singapore"
__license__ = "GPL2"
__credits__ = [""]

DEFAULT_SAMTOOLS_BINARY = "samtools"

# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


      
CoveragePeak = namedtuple('CoveragePeak', ['start', 'end', 'ori'])


def parse_peaks_file(peaks_file):
    """
    Parse a coverage peak file as produced by cov_peak_finder.py.
    Format is `start pos orientation`, where positions are zero-based,
    half-open and orientation is either 'F' or 'R'
    """

    peak_list = []
    fhandle = open(peaks_file, 'r')
    for line in fhandle:
        if line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue

        (start, end, ori) = line.split()
        start = int(start)
        end = int(end)
        
        assert start < end, ("Invalid position found in %s" % peaks_file)
        assert ori in ['F', 'R'], (
            "Invalid orientation '%s' found in %s" % (ori, peaks_file))

        peak = CoveragePeak(start=start, end=end, ori=ori)
        peak_list.append(peak)
    fhandle.close()
    
    return peak_list



def print_bam_header(bam_in, fh_out=sys.stdout, samtools='samtools'):
    """
    Print BAM header of `bam_in` to file handle `fh_out` using
    samtools binary `samtools`
    """

    cmd = '%s view -H %s' % (
        samtools, bam_in)

    #try:
    process = subprocess.Popen(cmd.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    #except OSError:
    #    import pdb; pdb.set_trace()
    (stdoutdata, stderrdata) =  process.communicate()

    retcode = process.returncode
    if retcode != 0:
        raise OSError("Called command exited with error code '%d'." \
                      " Command was '%s'. stderr was: '%s'" % (
                          retcode, cmd, stderrdata))
    for line in str.splitlines(stderrdata):
        if len(line.strip()) == 0:
            continue
        LOG.warn("samtools produced the following warning: %s" % line)

    for line in str.splitlines(stdoutdata):
        if len(line.strip()) == 0:
            continue
        fh_out.write(line + "\n")
    



def peak_remover(peaks_file, bam_in, fh_sam_out, samtools="samtools"):
    """
    FIXME
    """

    nremoved_fw = 0
    nremoved_rv = 0
    npassed = 0
    
    print_bam_header(bam_in, fh_sam_out, samtools)

    peaks = parse_peaks_file(peaks_file)

    rem_fw_start_pos = [p.start for p in peaks if p.ori == 'F']
    rem_rv_start_pos = [p.start for p in peaks if p.ori == 'R']

    LOG.info("Filtering forward reads starting at (zero-based): %s" % (rem_fw_start_pos))
    LOG.info("Filtering reverse reads starting at (zero-based): %s" % (rem_rv_start_pos))
    cmd = '%s view %s' % (
        samtools, bam_in)

    process = subprocess.Popen(cmd.split(),
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) =  process.communicate()

    retcode = process.returncode
    if retcode != 0:
        raise OSError("Called command exited with error code '%d'." \
                      " Command was '%s'. stderr was: '%s'" % (
                          retcode, cmd, stderrdata))

    for line in str.splitlines(stderrdata):
        if len(line.strip()) == 0:
            continue
        LOG.warn("samtools produced the following warning: %s" % line)
        
    for line in str.splitlines(stdoutdata):
        if len(line.strip()) == 0:
            continue
        line_split = line.split('\t')
        
        # SAM spec
        assert len(line_split) >= 11
        
        flag = int(line_split[1])
        pos = int(line_split[3])-1 # internally zero offset!
        
        if flag & 0x4 or flag & 0x200:
            LOG.warn("Detected unmapped fragment (0x4) or"
                     " not quality passed fragment: %s. Skippping" % (
                         line_split[0]))
            continue
        
        if flag & 0x10:
            # reverse
            if pos in rem_rv_start_pos:
                LOG.debug("Deleting (rv) line %s" % line)
                nremoved_rv += 1
                continue
        else:
            # forward
            if pos in rem_fw_start_pos:
                LOG.debug("Deleting (fw) line %s" % line)
                nremoved_fw += 1
                continue

        npassed += 1
        fh_sam_out.write(line + "\n")

    LOG.info("Num passed reads: %d. Forward/Reverse strand filter: %d/%d" % (
        npassed, nremoved_fw, nremoved_rv))
    # { fbam=CAGATC_1/CAGATC_1_remap_razers-i88.bam;
    #   samtools view -H $fbam;   samtools view $fbam | awk '{print $0}';
    # } | samtools view -bS - | samtools sort - schmock

    

def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" + \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    
    parser.add_option("-i", "--in",
                      dest="bam_in", # type="string|int|float",
                      help="Input BAM file)")
    parser.add_option("-o", "--out",
                      dest="sam_out", # type="string|int|float",
                      help="Output SAM file or '-' for stdout"
                       " (in which case you should piping through:"
                      " samtools view -bS - | samtools sort - outbase")
    parser.add_option("-p", "--peaks",
                      dest="peaks_file", # type="string|int|float",
                      help="File containig zero-based half open peak"
                      " description (as produced by cov_peak_finder.py)")
    parser.add_option("", "--samtools",
                      dest="samtools_binary", # type="string|int|float",
                      help="Optional: Path to samtools binary"
                      " (default: %s)" % DEFAULT_SAMTOOLS_BINARY)
    return parser



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

    
    if not opts.bam_in or not opts.sam_out or not opts.peaks_file:
        parser.error("Missing argument.")
        sys.exit(1)

    if not os.path.exists(opts.bam_in):
        LOG.fatal(
            "'%s' does not exist.\n" % opts.bam_in)
        sys.exit(1)

    if not os.path.exists(opts.peaks_file):
        LOG.fatal(
            "'%s' does not exist.\n" % opts.peaks_file)
        sys.exit(1)

    if opts.sam_out != "-" and os.path.exists(opts.sam_out):
        LOG.fatal(
            "Refusing to overwrite existing file '%s'.\n" % opts.sam_out)
        sys.exit(1)
        
    if not opts.samtools_binary:
        samtools_binary = DEFAULT_SAMTOOLS_BINARY
    else:
        samtools_binary = opts.samtools_binary

    if opts.sam_out == '-':
        fh_sam_out = sys.stdout
    else:
        fh_sam_out = open(opts.sam_out, 'w')
        
    peak_remover(opts.peaks_file,
                 opts.bam_in, fh_sam_out,
                 samtools_binary)

    if fh_sam_out != sys.stdout:
        fh_sam_out.close()

       
if __name__ == "__main__":
    main()
    LOG.info("Successful exit")
