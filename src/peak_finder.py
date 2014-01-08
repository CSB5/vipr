#!/usr/bin/env python
"""Finds and reports coverage peaks from PCR primers in a high
coverage mapping file containing reads coming from several PCR
fragments. The peaks are found by a simple three-sigma rule, which
works more or less. Reads found at these positions should be removed
before quality calibration and positions should also be ignored during
SNP calling.

See also peak_remover.py
"""


#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import gzip
import os
import datetime
import subprocess
import tempfile

#--- third-party imports
#
import numpy

#--- project specific imports
#
# /


MAX_LEN = 30000
MYNAME =  os.path.basename(sys.argv[0])
PLOT = False

timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

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



def median(values):
    """
    Taken from http://wiki.python.org/moin/SimplePrograms
    """
    
    copy = sorted(values)
    size = len(copy)
    assert size != 0, "Received empty list"

    if size % 2 == 1:
        retval = copy[(size - 1) / 2]
    else:
        retval = (copy[size/2 - 1] + copy[size/2]) / 2.0
    return float(retval)


    
def peaks_from_readstarts(fbam, samtools="samtools"):
    """
    FIXME
    """

    cmd = "%s view %s" % (samtools, fbam)
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
        LOG.warn("Got following stderr message from samtools: %s" % line)

    read_start_hist = dict()
    for ori in ['fw', 'rv']:
        read_start_hist[ori] = dict()

    readlen_hist = dict() # not strictly needed. only for backwards compatibilty
    
    for line in str.splitlines(stdoutdata):
        line_split = line.split("\t")

        # See SAM specification
        
        if len(line_split) < 11:
            continue

        flag = int(line_split[1])
        start_pos = int(line_split[3])-1 # internally zero offset!
        if flag & 0x4 or flag & 0x200:
            LOG.warn("Detected unmapped fragment (0x4) or"
                     " not quality passed fragment: %s. Skippping" % (
                         line_split[0]))
            continue

        if flag & 0x10:
            ori = 'rv'
        else:
            ori = 'fw'

        # increment dict value if it exists, otherwise set to 1
        hist = read_start_hist[ori]
        hist[start_pos] = hist.get(start_pos, 0) + 1

        # not strictly needed (only plotting). heuristic anyway
        cigar = line_split[5]
        if cigar.endswith('M'):
            try:
                readlen = int(cigar[:-1])
                readlen_hist[readlen] = readlen_hist.get(readlen, 0) + 1
            except ValueError:
                pass
        
    peak_pos = dict()
    for ori in ['fw', 'rv']:
        mean = numpy.mean(read_start_hist[ori].values())
        std = numpy.std(read_start_hist[ori].values())
        thresh = mean + 3*std
        peak_pos[ori] = [k for (k, v) in
                           read_start_hist[ori].iteritems() if v>thresh]

        
    # paranoia (and again: only for backward compatibility)
    if len(readlen_hist.values()) != 1:
        LOG.warn("Detected reads with different length!")
    readlen = sorted(readlen_hist.items(), key=lambda x: x[1])[-1][0]

    return ([(s, s+readlen) for s in peak_pos['fw']],
            [(s, s+readlen) for s in peak_pos['rv']],
            readlen)



def genome_from_bam(fbam, fh_genome, samtools="samtools"):
    """
    Create A bedtools 'genome' file

    Arguments:
    - `fbam`: bam file to create genome file from
    - `fh_genome`: bedtools genome file handle
    - `samtool`: path to samtools binary
    """
    
    cmd = "%s view -H %s" % (samtools, fbam)
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
        LOG.warn("Got following stderr message from samtools: %s" % line)

    # local requirement would be == 1
    #if stdoutdata.count("@SQ") != 1:
    #     raise ValueError(
    #         "Did not find exactly one @SQ line in BAM file '%s'"  % (
    #             fbam))
     
    for line in str.splitlines(stdoutdata):
        if line.startswith("@SQ"):
            line_split = line.split('\t')
            chrname = line_split[1]
            assert chrname.startswith("SN:")
            chrname = chrname[3:]

            chrlen = line_split[2]
            assert chrlen.startswith("LN:")
            chrlen = int(chrlen[3:])
            
            fh_genome.write("%s\t%d\n" % (chrname, chrlen))

    
def parse_coverage_bam(fbam, genomecoveragebed="genomeCoverageBed"):
    """
    
    Arguments:
    - `fbam`: file to read coverage from. Needs samtools and bedtools (genomeCoverageBed)
    - `genomecoveragebed`: path to genomeCoverageBed binary
    """

  
    # forward and reverse coverage
    # int is essential; see note below
    fw_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)
    rv_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)

    
    file_genome = tempfile.NamedTemporaryFile(delete=False)
    genome_from_bam(fbam, file_genome)
    file_genome.close()

    basic_cmd = "%s -ibam %s -g %s -d" % (
        genomecoveragebed, fbam, file_genome.name)



    for strand in ["+", "-"]:
        cmd = "%s -strand %s" % (basic_cmd, strand)
        try:
            process = subprocess.Popen(cmd.split(),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        except OSError:
            raise OSError, ("It seems %s is not installed" % cmd.split()[0])
        (stdoutdata, stderrdata) =  process.communicate()
        retcode = process.returncode
        if retcode != 0:
            raise OSError("Called command exited with error code '%d'." \
                          " Command was '%s'. stderr was: '%s'" % (
                              retcode, cmd, stderrdata))
        for line in str.splitlines(stderrdata):
            if len(line.strip()) == 0:
                continue
            LOG.warn("Got following stderr message from genomeCoverageBed: %s" % line)
        for line in str.splitlines(stdoutdata):
            if len(line) == 0:
                continue
            (chr, pos, cov) = line.split('\t')
            pos = int(pos)-1 # we use zero offset
            cov = int(cov)

            assert pos < MAX_LEN

            if strand == '+':
                fw_cov[pos] = cov
            elif strand == '-':
                rv_cov[pos] = cov

    os.unlink(file_genome.name)

    return (numpy.trim_zeros(fw_cov, trim='b'),
            numpy.trim_zeros(rv_cov, trim='b'))


    
def parse_coverage_razers(fhandle):
    """
    
    Arguments:
    - `fhandle`: file handle to read razers mapping from
    """

    # forward and reverse coverage
    # int is essential; see note below
    fw_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)
    rv_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)

    readlen = None
    
    for line in fhandle:
        if line.startswith('#'):
            continue     
        #print "DEBUG %s" % line,
        
        # taken from razers.py
        line_split = line.split('\t')
        assert len(line_split) == 8 or len(line_split) == 11, (
            "Was expecting 8 or 11 elements but got %d in line: %s" % (len(line_split), line_split))
        #if len(line_split) > 8:
        #    raise ValueError, (
        #        "Not tested with paired end reads")
        (rname,
         rbegin,
         rend,
         gstrand,
         gname,
         gbegin,
         gend,
         percid) = line_split[:8]
        # don't need pairing information
        
        # razers position format: zero-based & half open. as python
        # slices/ranges
        gbegin = int(gbegin)
        gend = int(gend)
        
        rbegin = int(rbegin)
        rend = int(rend)
        if readlen:
            assert rend-rbegin == readlen, (
                "Can't handle reads of different lengths")
        else:
            readlen = rend-rbegin
            
        # this will work even with indels. however, the indel position
        # will also be treated as covered! to prevent this, we'd have
        # to parse the sequences as well which is not implemented at
        # the moment
        assert gend-gbegin == readlen, (
            "Can't handle indels")

        if gstrand == 'F':
            fw_cov[gbegin:gend] += 1
        elif gstrand == 'R':
            rv_cov[gbegin:gend] += 1
        else:
            raise ValueError, "Unknown strand orientation '%s'" % gstrand

    return (numpy.trim_zeros(fw_cov, trim='b'),
            numpy.trim_zeros(rv_cov, trim='b'),
            readlen)



def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog:" + __doc__ +  \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    parser.add_option("-b", "--bam",
                      dest="fbam", # type="string|int|float"
                       help="BAM input file (needs samtools installed)")
    parser.add_option("-o", "--output",
                      dest="fout", # type="string|int|float"
                       help="peak file (output)")
    parser.add_option("-i", "--razers",
                      dest="frazers", # type="string|int|float"
                       help="Razers input file (instead of BAM. Can be zipped; use '-' for stdin)")
    parser.add_option("-p", "--plot",
                      dest="fplot", # type="string|int|float"
                       help="Optional: plot file for visualising found positions (output)")
    #parser.add_option("-r", "--readlen",
    #                  dest="readlen", type="int",
    #                   help="Read length (needed for BAM input)")
    return parser



def main():
    """
    The main function
    """


    
    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    #if len(args):
    #    parser.error("Unrecognized arguments found: %s." % (
    #        ' '.join(args)))
    #    sys.exit(1)
        
    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    close_handle = True

    if opts.fbam and opts.frazers:
        parser.error("Can only use either RazerS or BAM as input.")
        sys.exit(1)
        
    if not opts.fbam and not opts.frazers:
        parser.error("input file argument missing.")
        sys.exit(1)

    if opts.frazers and opts.frazers != "-" and not os.path.exists(opts.frazers):
        LOG.fatal(
            "input file '%s' does not exist.\n" % opts.frazers)
        sys.exit(1)

    if opts.fbam and not os.path.exists(opts.fbam):
        LOG.fatal(
            "input file '%s' does not exist.\n" % opts.fbam)
        sys.exit(1)
    #if opts.fbam and not opts.readlen:
    #    LOG.fatal(
    #        "Need readlen when using a BAM files as input. Try 'samtools view %s | awk \'{readlen[length($10)]+=1} END {for (word in readlen) {print word}}\'.\n" % opts.fbam)
    #    sys.exit(1)
        
        
    if not opts.fout:
        parser.error("output file argument missing.")
        sys.exit(1)
    if os.path.exists(opts.fout):
        LOG.fatal(
            "Refusing to overwrite already existing output file '%s'.\n" % opts.fout)
        sys.exit(1)

    if opts.frazers:
        if opts.frazers == "-":
            fhandle = sys.stdin
            close_handle = False
        elif opts.frazers[-3:] == ".gz":
            fhandle = gzip.open(opts.frazers, 'r')
        else:
            fhandle = open(opts.frazers, 'r')

        LOG.warn("This uses legacy code. Use BAM input instead")
        (fw_cov, rv_cov, readlen) = parse_coverage_razers(fhandle)
     
        if close_handle:
            fhandle.close()

        # make zeros ones to prevent division by zero errors
        fw_cov[numpy.where(fw_cov==0)] = 1
        rv_cov[numpy.where(rv_cov==0)] = 1
     
        # fw_*[pos] has the change coming from forward and
        # rv_*[pos] has the change coming from reverse.
        # holds rv_cov[pos] - rv_cov[pos+1]
        #
        # overhanging ends are treated as zeros/ones.
     
        # fc's need to be  float
        fw_fc = fw_cov.astype(numpy.float) / numpy.append(
            numpy.ones((1)), fw_cov[:-1].astype(numpy.float))
     
        rv_fc = rv_cov.astype(numpy.float) / numpy.append(
            numpy.ones((1)), rv_cov[:-1].astype(numpy.float))
        
        dev = 0.5
        fw_peaks_up = numpy.where(fw_fc > 1 + dev)[0]
        fw_peaks_down = numpy.where(fw_fc < 1 - dev)[0]
     
        rv_peaks_up = numpy.where(rv_fc > 1 + dev)[0]
        rv_peaks_down = numpy.where(rv_fc < 1 - dev)[0]
     
        primer_peaks_fw = []
        for p_up in fw_peaks_up:
            for p_down in fw_peaks_down:
                if abs(p_down - p_up) == readlen:
                    primer_peaks_fw.append((p_up, p_down))
     
        primer_peaks_rv = []
        for p_up in rv_peaks_up:
            for p_down in rv_peaks_down:
                if abs(p_down - p_up) == readlen:
                    primer_peaks_rv.append((p_up, p_down))
                    
        #if len(primer_peaks_fw) != len(primer_peaks_rv):
        #    LOG.warn("Number of forward and reverse peaks does not match (seem peaks undetected?)")


    elif opts.fbam:
        # coverage only needed for plotting
        (fw_cov, rv_cov) = parse_coverage_bam(opts.fbam)
        (primer_peaks_fw, primer_peaks_rv,
         readlen) = peaks_from_readstarts(opts.fbam)


    else:
        LOG.fatal("Internal error: got neither BAM nor RazerS file")
        sys.exit(1)


         
    #import pdb; pdb.set_trace()
    
    fhandle = open(opts.fout, 'w')
    fhandle.write("# Created with %s at %s\n" % (
        MYNAME, timestamp))
    fhandle.write("# Detected peaks with plateau of length %d:\n" % readlen)
    fhandle.write("# Reported as: <start> <end> <orientation>\n")
    fhandle.write("# where start-end range is zero-based, half-open and orientation=[FR]\n")
    for p in primer_peaks_fw:
        fhandle.write("%d %d %s\n" % (p[0], p[1], 'F'))
    for p in primer_peaks_rv:
        fhandle.write("%d %d %s\n" % (p[0], p[1], 'R'))
    fhandle.close()
                      
    #import pdb; pdb.set_trace()

    if opts.fplot:
        import matplotlib.pyplot as pyplot
        fig = pyplot.figure()
        ax1 = pyplot.subplot(111)
        
        #ax1.set_yscale('log');# breaks broken_barh

        # plot coverage
        #
        ax1.plot(fw_cov, color='green')
        ax1.plot(rv_cov, color='red')

        # plot peaks
        #
        (bbar_ymin, bbar_height) = (0, max(max(fw_cov), max(rv_cov))/10.0)
        ax1.broken_barh([(p[0], p[1]-p[0]) for p in primer_peaks_fw],
                        (bbar_ymin, bbar_height), color="blue")
        ax1.broken_barh([(p[0], p[1]-p[0]) for p in primer_peaks_rv],
                        (bbar_ymin, bbar_height), color="orange")

        count = 0       
        for p in primer_peaks_fw:
            count += 1
            ax1.annotate('%d' % p[0],
                         xy = (p[0], bbar_height*(1.25+(0.25 * (count%2)))),
                         horizontalalignment='left',
                         verticalalignment='center',
                         fontsize=6, color="blue")

        for p in primer_peaks_rv:
            count += 1
            ax1.annotate('%d' % p[0],
                         xy = (p[0], bbar_height*(1.25+(0.25 * (count%2)))),
                         horizontalalignment='right',
                         verticalalignment='center',
                         fontsize=6, color="orange")
                         #xy = (p[0] + (p[1]-p[0])/2.0, bbar_height),
                         #xytext = (p[0] + (p[1]-p[0])/2.0 + max(fw_cov)/25, bbar_height*(1.5+count%1)),
                         #arrowprops=dict(facecolor='black', shrink=0.05),
        
        # undetected forward/reverse in read-GATCAG__ref-05K2928DK1-Den1__razers-rr100-i85-unique.out.gz
        #
        #ax1.broken_barh([(3159-1, READLEN), (3416-1, READLEN)],
        #                (bbar_ymin, bbar_height/2), facecolors="brown")
        #ax1.broken_barh([(3417-26, READLEN), (3761-26, READLEN)],
        #                (bbar_ymin, bbar_height/2), facecolors="yellow")

        #ax1.axis([1000, 5000, 0, 5000])

        # write to file
        #
        fileext = os.path.splitext(opts.fplot)[1]
        pyplot.savefig(opts.fplot, format=fileext[1:])
        #sys.stderr.write("Pseudointeractivity...\n")
        #import pdb; pdb.set_trace();

    
if __name__ == "__main__":
    main()
    LOG.info("successful exit")
