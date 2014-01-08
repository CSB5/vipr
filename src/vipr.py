#!/usr/bin/env python
"""VIral Pipeline Runner
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
#import tempfile
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import subprocess
import datetime

#--- third-party imports
#
#/

#--- project specific imports
#
from Bio import SeqIO


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013 Genome Institute of Singapore"
__license__ = "GPL2"
__credits__ = [""]
__status__ = "eternal alpha"


DEFAULT_NUM_THREADS = 2


# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


class ViprMakefile(object):
    """Interface for makefile creation for Viral pipeline
    """

    myname = os.path.basename(sys.argv[0])

    
    def __init__(self, makefile_name, outprefix,
                 read_file_list, illumina_enc, 
                 ref_fa, primer_fa, 
                 snv_pred_against_self=1,
                 primer_len=25, num_threads=2):
        """init function
        """

        self.makefile_name = makefile_name
        assert not os.path.exists(makefile_name)
        self.makefile_fh = None
        self.outprefix = outprefix

        # all filename need to be relative to filename
        self.read_file_list = []
        for f in read_file_list:
            self.read_file_list.append(
                os.path.relpath(f, os.path.dirname(makefile_name)))
        self.ref_fa = os.path.relpath(ref_fa, 
                                      os.path.dirname(makefile_name))
        self.primer_fa = os.path.relpath(primer_fa,
                                         os.path.dirname(makefile_name))
        self.snv_pred_against_self = snv_pred_against_self
        self.primer_len = primer_len
        self.num_threads = num_threads
        self.illumina_enc = illumina_enc



    def write_header(self):
        """Writes makefile header consisting only of comments
        """
        
        timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        self.makefile_fh.write("# Automatically created by %s on %s\n"
        "#\n"
        "# Makefile resources:\n"
        "# - http://en.wikipedia.org/wiki/Makefile\n"
        "# - http://www.cprogramming.com/tutorial/makefiles.html\n"
        "# Makefiles in Bioinfo:\n"
        "# - http://www.bioinformaticszen.com/post/decomplected-workflows-makefiles/\n"
        "# - http://archive.nodalpoint.org/2007/03/18/a_pipeline_is_a_makefile\n"
        "# - http://janos.binder.hu/?p=190\n"
        "#\n"
        "# Advanced Makefile Tricks\n"
        "# - http://www.cprogramming.com/tutorial/makefiles_continued.html\n"
        "# Makefile notes\n"
        "# - http://www.acsu.buffalo.edu/~charngda/makefile.html\n"
        "\n" % (self.myname, timestamp))


    def write_userdata(self):
        """Writes the actual user-defined variables to Makefile"""

        self.makefile_fh.write("# Predefined data\n")
        self.makefile_fh.write("#\n")
        self.makefile_fh.write("PREFIX := %s\n" % self.outprefix)
        self.makefile_fh.write("S1 := %s\n" % self.read_file_list[0])
        self.makefile_fh.write("# for PE define S2; for SR leave S2 blank\n")
        if len(self.read_file_list)==2:
            self.makefile_fh.write("S2 := %s\n" % self.read_file_list[1])
        else:
            self.makefile_fh.write("S2 := \n")
        self.makefile_fh.write("INIT_REF := %s\n" % self.ref_fa)
        self.makefile_fh.write("SNV_PRED_AGAINST_SELF := %d\n" % (
            self.snv_pred_against_self))

        self.makefile_fh.write("PRIMER_FILE := %s\n" % self.primer_fa)
        self.makefile_fh.write("NUM_THREADS := %d\n" % self.num_threads)
        self.makefile_fh.write(
            "ILLUMINA := %d\n" % (1 if self.illumina_enc else 0))
        self.makefile_fh.write("PRIMER_LEN := %d\n" % self.primer_len)
        self.makefile_fh.write("\n")

        
    def write_spacer(self, num=1):
        """Add spacer to Makefile for readability"""
        for i in range(num):
            self.makefile_fh.write("\n")


    def write_predef(self):
        """Writes the predefined bit i.e. the actual workflow to the
         Makefile"""

        predef = """MAP_REF_FA := $(PREFIX)_mapping_ref.fa
CONS_FA := $(PREFIX)_cons.fa
CONS_FA_MASKED := $(CONS_FA:%.fa=%_masked.fa)
PRIMER_POS_CONS_FA := $(CONS_FA:%.fa=%_primer-pos.txt)
BWA_UNIQ_BAM := $(PREFIX)_bwa-uniq.bam
MAPPING_SUCCESS := $(BWA_UNIQ_BAM:%.bam=%-mapping-success.txt)
MAPPING_COVPLOT := $(BWA_UNIQ_BAM:%.bam=%-coverage.pdf)
PRIMER_POS := $(MAP_REF_FA:%.fa=%_primer-pos.txt)
INCLUDE_BED := $(MAP_REF_FA:%.fa=%_include.bed)
DUPS_MARKED_BAM := $(BWA_UNIQ_BAM:%.bam=%-dups-marked.bam)
RECAL_BAM := $(DUPS_MARKED_BAM:%.bam=%-recal.bam)
LOFREQ_RAW := $(RECAL_BAM:%.bam=%_lofreq-raw.snp)
LOFREQ_FINAL := $(LOFREQ_RAW:%raw.snp=%final.snp)
SNP_CLUSTER := $(LOFREQ_FINAL:%.snp=%_cluster.txt)

ILL_ENC_ARG = 
ifeq ($(ILLUMINA),1)
\tILL_ENC_ARG = --illumina
endif
ifdef S2
\tREADS_ARG = -f $(S1) -g $(S2)
\tMAPPING_SUCCESS_PE_ARG = -p
else
\tREADS_ARG = -f $(S1)
\tMAPPING_SUCCESS_PE_ARG =
endif

$(CONS_FA): $(S1) $(S2) $(INIT_REF) $(PRIMER_FILE)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tbam2cons_iter.sh $(ILL_ENC_ARG) $(READS_ARG) -r $(INIT_REF) -t $(NUM_THREADS) --force -o .$@ > $@.log 2>&1
\tprimer_pos_from_seq.sh -p $(PRIMER_FILE) -r .$@ --force -o $(PRIMER_POS_CONS_FA) >> $@.log 2>&1
\tmask_primer.py --force -i .$@ -o $(CONS_FA_MASKED) -p $(PRIMER_POS_CONS_FA) >> $@.log 2>&1
\tmv .$@ $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(MAP_REF_FA): $(S1) $(S2) $(INIT_REF) $(CONS_FA)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\t@if [ $(SNV_PRED_AGAINST_SELF) -eq 1 ]; then\\
\t  cp $(CONS_FA) $@;\\
\telse\\
\t  cp $(INIT_REF) $@;\\
\tfi
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(BWA_UNIQ_BAM): $(S1) $(S2) $(CONS_FA)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tbwa_unique.sh $(ILL_ENC_ARG) $(READS_ARG) -r $(MAP_REF_FA) -t $(NUM_THREADS) -o .$@  > $@.log 2>&1
\tmv .$@ $@
\tsamtools index $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(MAPPING_SUCCESS): $(S1) $(BWA_UNIQ_BAM)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tcoverage_plot.py --force --log $(MAPPING_COVPLOT).txt -o $(MAPPING_COVPLOT) -b $(BWA_UNIQ_BAM)
\tmapping_success.sh $(MAPPING_SUCCESS_PE_ARG) -f $(S1) -b $(BWA_UNIQ_BAM) > .$@
\tmv .$@ $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

# FIXME if PRIMER_FILE is missing use peak-caller primer_pos_from_peaks.py
$(PRIMER_POS): $(PRIMER_FILE) $(MAP_REF_FA)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tprimer_pos_from_seq.sh -p $(PRIMER_FILE) -r $(MAP_REF_FA) --force -o .$@ > $@.log 2>&1
\tmv .$@ $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(DUPS_MARKED_BAM): $(BWA_UNIQ_BAM) $(PRIMER_POS)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tmark_primer.py --primer-len $(PRIMER_LEN) -i $(BWA_UNIQ_BAM)  -p $(PRIMER_POS) --force -o .$@ > $@.log 2>&1
\tmv .$@ $@
\tsamtools index $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(RECAL_BAM): $(DUPS_MARKED_BAM) $(MAP_REF_FA)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tbase_qual_calib_wrapper.sh -t $(NUM_THREADS) -i $(DUPS_MARKED_BAM) -r $(MAP_REF_FA) -o .$@  > $@.log 2>&1
\tmv .$@ $@
\tsamtools index $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(LOFREQ_FINAL): $(PRIMER_POS) $(RECAL_BAM) $(MAP_REF_FA)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tprimer_pos_to_bed.py --force -i $(PRIMER_POS) -p $(PRIMER_LEN) -b $(RECAL_BAM) -o $(INCLUDE_BED) > $@.log 2>&1
\tlofreq_snpcaller.py --force -f $(MAP_REF_FA) -b $(RECAL_BAM) -o $(LOFREQ_RAW) -l $(INCLUDE_BED)  >> $@.log 2>&1
\tlofreq_filter.py --force --strandbias-holmbonf --min-cov 10 -i $(LOFREQ_RAW) -o .$@ >> $@.log 2>&1
\tmv .$@ $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"

$(SNP_CLUSTER): $(LOFREQ_FINAL)
\t@echo "INFO: Starting to work on $@ (`date "+%Y-%m-%d %H:%M"`)"
\t@echo "INFO: Check $@.log if things fail."
\tlofreq2_cluster.py -i $(LOFREQ_FINAL) -o .$@ > $@.log 2>&1
\tmv .$@ $@
\t@echo "INFO: Done (`date "+%Y-%m-%d %H:%M"`)"


all: $(CONS_FA) $(MAP_REF_FA) $(MAPPING_SUCCESS) $(LOFREQ_FINAL) $(SNP_CLUSTER)
        """

        self.makefile_fh.write(predef)


            
    def write(self):
        """Writes the Makefile
        """

        self.makefile_fh = open(self.makefile_name, "w")
        self.write_header()
        self.write_spacer()
        self.write_userdata()
        self.write_spacer()
        self.write_predef()
        self.write_spacer()
        self.makefile_fh.close()




def which(prog):
    """make sure prog can be run
    """

    try:
        subprocess.call([prog],
                        stderr=subprocess.PIPE,
                        stdout=subprocess.PIPE,)
        return True
    except OSError:
        return False



def run_fastqc(fastq_file_list, out_dir=None,
               num_threads=DEFAULT_NUM_THREADS):
    """Wrapper for running fastqc"""

    assert isinstance(fastq_file_list, list), (
        "Need fastq files as list")

    all_fastqc_files_exist = True
    for fastq in fastq_file_list:
        if not os.path.exists(
                derive_fastqc_data_filename(fastq, out_dir)):
            all_fastqc_files_exist = False
    if all_fastqc_files_exist:
        LOG.info("Reusing fastqc files")
        return

    cmd_list = ['fastqc']
    cmd_list.append("-q")
    cmd_list.extend(["-t", "%d" % num_threads])
    if out_dir:
        cmd_list.extend(["-o", "%s" % out_dir])
    cmd_list.extend(fastq_file_list)
    try:
        LOG.info("Executing: %s" % ' '.join(cmd_list))
        subprocess.check_call(cmd_list)
    except OSError:
        LOG.fatal("Couldn't execute: %s" % (
            ' '.join(cmd_list)))
        raise
    except subprocess.CalledProcessError:
        LOG.fatal("%s failed: %s" % (
            cmd_list[0], ' '.join(cmd_list)))
        raise

    
def derive_fastqc_data_filename(fastq, outdir=None):
    """Derive the FastQC data file from fastq input on which FastQC
    was run. Set outdir if fastqc was run with -o. No checks will be
    done."""

    # remove gz first and then remaining extension
    base = os.path.splitext(fastq.rsplit(".gz")[0])[0]
    if outdir:
        base = os.path.join(outdir, os.path.basename(base))
    return os.path.join(base + "_fastqc", "fastqc_data.txt")



def enc_from_fastqc_data(fastqc_data_filename):
    """Derive encodig from FastQC data file

    Return values are as defined in
    http://www.biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-module.html:

    "qual" means simple quality files using PHRED scores (e.g. from
    Roche 454)
    
    "fastq" means Sanger style FASTQ files using PHRED scores and an
    ASCII offset of 33 (e.g. from the NCBI Short Read Archive and
    Illumina 1.8+). These can potentially hold PHRED scores from 0 to
    93.
    
    "fastq-sanger" is an alias for "fastq".
    
    "fastq-solexa" means old Solexa (and also very early Illumina)
    style FASTQ files, using Solexa scores with an ASCII offset 64.
    These can hold Solexa scores from -5 to 62.

    "fastq-illumina" means newer Illumina 1.3 to 1.7 style FASTQ
    files, using PHRED scores but with an ASCII offset 64, allowing
    PHRED scores from 0 to 62.
    """

    fh = open(fastqc_data_filename, 'r')
    for line in fh:
        if line.startswith("Encoding"):
            enc_str = line.split("\t")[1]
            break
    fh.close()

    if not enc_str:
        raise ValueError, (
            "Couldn't find encoding in %s" % (fastqc_data_filename))
    if "Illumina 1.5" in enc_str:
        return "fastq-illumina"
    elif "Sanger" in enc_str:
        return "fastq-sanger"
    else:
        raise ValueError, (
            "Don't know about encoding called '%s'" % enc_str)

    
def enc_from_fastq(fastq, outdir=None, num_threads=1):
    """Will derive quality encoding from fastq. Will reuse fastqc info
    in fastq-dir or outdir if existant and run fastqc in outdir
    otherwise
    """

    if os.path.exists(derive_fastqc_data_filename(fastq)):
        enc = enc_from_fastqc_data(
            derive_fastqc_data_filename(fastq))
    elif os.path.exists(derive_fastqc_data_filename(fastq, outdir)):
        enc = enc_from_fastqc_data(
            derive_fastqc_data_filename(fastq, outdir))
    else:
        run_fastqc([fastq], outdir, num_threads)
        enc = enc_from_fastqc_data(
            derive_fastqc_data_filename(fastq, outdir))
    return enc


def primer_len_from_fasta(primer_fa):
    """Derive maximum primer length from fasta file"""

    try:
        seq_recs = list(SeqIO.parse(open(primer_fa, 'r'), 'fasta'))
    except:
        LOG.fatal("Doesn't look like %s is a valid fasta file" % (primer_fa))
        sys.exit(1)

    assert len(seq_recs)>0, (
        "No sequences found in %s" % (primer_fa))
    max_len = max([len(sr.seq) for sr in seq_recs])
    assert max_len > 0, (
        "Primer length is 0!? Are you sure %s is a valid fasta file?" % primer_fa)
    assert max_len < 100, (
        "Primer length is > 100!?")

    return max_len


def drop_readme(readme_filename, result_prefix):
    """Add a README to output directory explaining files"""

    readme_fh = open(readme_filename, 'w')

    # note: keep this in sync with the actual produced files
    
    template = """This file lists the most important result files in this directory.

Consensus/Dominant/Master sequence
(NOTE: primer regions might differ between samples. If in doubt use the masked version; see below)
__PREFIX___cons.fa

Consensus/Dominant/Master sequence with masked primer regions:
__PREFIX___cons_masked.fa

Sequence used as reference for mapping / low-freq. SNV prediction:
__PREFIX___mapping_ref.fa

Primer positions detected in mapping reference:
__PREFIX___mapping_ref_primer-pos.txt

Raw, uncalibrated mapping (BWA; uniquely mapped reads):
__PREFIX___bwa-uniq.bam

Mapping stats of raw mapping (__PREFIX___bwa-uniq.bam):
__PREFIX___bwa-uniq-mapping-success.txt

Coverage plot for raw mapping (__PREFIX___bwa-uniq.bam):
__PREFIX___bwa-uniq-coverage.pdf

Final mapping, with primer dups marked/clipped and qualities recalibrated:
__PREFIX___bwa-uniq-dups-marked-recal.bam

Regions analyzed by LoFreq (excluding primers):
__PREFIX___mapping_ref_include.bed

Filtered, low-frequency variants (for e.g. csv-import to Excel):
__PREFIX___bwa-uniq-dups-marked-recal_lofreq-final.snp
"""

    readme_fh.write(template.replace("__PREFIX__", result_prefix))

    readme_fh.close()


def call_make():
    """Calls make and redirects stdout and stderr to logfile. Will
    raise exception if make failed
    """

    LOG.info("Running in make in %s. This will take a while.." % os.getcwd())
    log_filename = "Makefile.log"

    log = open(log_filename, 'w')
    try:
        subprocess.check_call(['make', 'all'], 
                              stdout=log, stderr=log)
    except:
        LOG.fatal("Analysis in %s failed. Check %s,"
                  " fix error and run again with"
                  " --continue" % (os.getcwd(), log_filename))
        sys.exit(1)
    log.close()

    # the following was supposed to mimic tee but doesn't work at all
    
#    proc = subprocess.Popen(["make", "all"], 
#                            stdout=subprocess.PIPE, 
#                            stderr=subprocess.PIPE)
#    with open(log_filename, 'w') as log:
#        while proc.poll() is None:
#            line = proc.stderr.readline()
#            if line:
#                sys.stderr.write(line)
#                log.write(line)
#            line = proc.stdout.readline()
#            if line:
#                sys.stdout.write(line)
#                log.write(line)


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "usage: %prog [options]\nType %prog -h to get more help"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="Optional: be verbose")
    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="Optional: enable debugging")
    #parser.add_option("-l", "--lib-dir",
    #                  dest="libdir",
    #                  help="Library directory containing fastq-files."
    #                  " You can give multiple dirs (keep in quotes)"
    #                  " Example: /mnt/SolexaPool/HS002/HS002-PE-R00053_BC0JPFACXX/output_HS002-PE-R00053_CASAVA-v1.8.2_20120810/CASAVA-v1.8.2_Unaligned_20120810/Project_MUX216")
    #parser.add_option("-m", "--mux-dir",
    #                  dest="muxdir",
    #                  help="MUX directory containing library sub-directories."
    #                  " If given, will analyse all library sub-directories."
    #                  " Example: /mnt/SolexaPool/HS002/HS002-PE-R00053_BC0JPFACXX/output_HS002-PE-R00053_CASAVA-v1.8.2_20120810/CASAVA-v1.8.2_Unaligned_20120810/Project_MUX216/Sample_WEB016")
    #parser.add_option("-o", "--out-dir",
    #                  dest="outdir_root",
    #                 default=OUT_DIR_ROOT,
    #                  help="Output directory to dump results in"
    #                  " (defaults to %s followed by mux or"
    #                  " sample-name)" % (OUT_DIR_ROOT))
    parser.add_option("-r", "--ref-fa",
                      dest="ref_fa",
                      help="Reference fasta file")
    parser.add_option("-1", "--s1",
                      dest="s1",
                      help="1st sequence read file.")
    parser.add_option("-2", "--s2",
                      dest="s2",
                      help="2nd sequence read file, if paired-end")
    parser.add_option("-p", "--primers",
                      dest="primer_fa",
                      help="Fasta file containing primer sequences")
    parser.add_option("-o", "--outdir",
                      dest="out_dir",
                      help="Output directory")
    parser.add_option("", "--outprefix",
                      dest="outprefix",
                      help="Prefix for all output filenames"
                      " (defaults to basename of read-file)")
    parser.add_option("", "--snvs-against-ref",
                      dest="snvs_against_ref",
                      action="store_true",
                      help="Predict low-freq. SNVs against given ref,"
                      " not the sample consensus")
    parser.add_option("-t", "--threads",
                      dest="num_threads",
                      type="int",
                      default=DEFAULT_NUM_THREADS,
                      help="Whenvever possible use this number of"
                      " processors (default %d)" % DEFAULT_NUM_THREADS)
    parser.add_option("-c", "--continue",
                      dest="continue_interrupted",
                      action="store_true",
                      help="Continue interrupted analysis in outdir."
                      " Needs only outdir as other argument")
    parser.add_option("", "--exit-before-make",
                      dest="exit_before_make",
                      action="store_true",
                      help="EXPERT ONLY: Exit before executing make")

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

    #if opts.verbose:
    #    LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    # file check
    for (filename, descr, direction, mandatory) in [
            (opts.s1, "Sequence reads (1st file)", 'in', True),
            (opts.s2, "Sequence reads (2nd file)", 'in', False),
            (opts.primer_fa, "Primer fasta file", 'in', True),
            (opts.ref_fa, "Reference fasta file", 'in', True),
            ]:

        # exception:
        # none of the input files are needed if --continue was given
        if opts.continue_interrupted:
            break
        
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

    # logic checks of arguments
    #
    #if not opts.continue_interrupted and os.path.exists(opts.out_dir):
    #    LOG.fatal("Directory %s already exists" % (opts.out_dir))
    #    sys.exit(1)
    if not opts.out_dir:
        LOG.fatal("Missing output directory argument")
        sys.exit(1)
    if opts.continue_interrupted:
        if not os.path.exists(opts.out_dir):
            LOG.fatal("Can't continue in directory %s because it"
                      " doesn't exist" % (opts.out_dir))
            sys.exit(1)

        if opts.exit_before_make:
            LOG.info("exit before make (to continue manually: cd %s && make all > Makefile.log 2>&1)" % os.getcwd())
            return
        os.chdir(opts.out_dir)
        call_make()
        return

    if opts.s2:
        assert opts.s1 != opts.s2

    if not opts.outprefix:
        outprefix = os.path.basename(os.path.splitext(
            opts.s1.rsplit(".gz")[0])[0])
    else:
        outprefix = opts.outprefix
        
    primer_len = primer_len_from_fasta(opts.primer_fa)

    if not os.path.exists(opts.out_dir):
        os.mkdir(opts.out_dir)
    
    read_file_list = [opts.s1]
    if opts.s2:
        read_file_list.append(opts.s2)

    enc = enc_from_fastq(opts.s1, opts.out_dir, opts.num_threads)
    illumina_enc = False
    if 'solexa' in enc or 'illumina' in enc:
        illumina_enc = True

    makefile_name =  os.path.join(opts.out_dir, "Makefile")
    if os.path.exists(makefile_name):
        LOG.fatal("I'm confused. Makefile %s already exists"
                  " but continue option was not given."
                  " If I should just continue call me with --continue."
                  " If you want me to start anew, delete just %s"
                  " (or the whole output directory)" % (
                      makefile_name, makefile_name))
        sys.exit(0)

    makefile_obj = ViprMakefile(
        makefile_name, outprefix, read_file_list, 
        illumina_enc, opts.ref_fa, opts.primer_fa, 
        not opts.snvs_against_ref, primer_len, opts.num_threads)    
    makefile_obj.write()

    readme_filename = os.path.join(opts.out_dir, "README.txt")
    drop_readme(readme_filename, outprefix)

    if opts.exit_before_make:
        LOG.info("exit before make (to continue manually: cd %s && make all > Makefile.log 2>&1)" % os.getcwd())
        return
    os.chdir(opts.out_dir)
    call_make()


if __name__ == "__main__":
    main()
    LOG.info("Successful exit")
    # FIXME normalize path for INIT_REF and PRIMER_FILE
    # use os.path.abspath but careful after cwd()


