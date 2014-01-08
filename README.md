ViPR
====

This is a collection of scripts used for the analysis of RNA viruses,
based on LoFreq. The main wrapper is vipr.py ("Viper" written in
Python though :) which will create a makefile defining the pipeline
setps and will subsequently call make. Some of these scripts are quite
specific to small genomes, others you might able to reuse some of them
for larger - e.g. bacterial - genomes.

*WARNING: Use at your own risk*: this pipeline is not based on a
proper pipeline framework (but is being ported to ruffus at the
moment), uses some outdated scripts, has a huge dependency list and
might make assumptions that are not met in your setup.

Installation
------------

The impatient can simply use the full path for calling the scripts in
./src or add that directory to their PATH variable.

For a proper
installation, cd into this directory and then use

    $ python setup.py
    
install --prefix for installing the scripts.

For more options see

    $ python setup.py install -h

Note: if you're using a special prefix for installation, please make
sure that the corresponding bin directory is in your PATH and the
corresponding python library path is in your PYTHONPATH.


Dependencies
------------

A (very likely) incomplete list of dependencies:
- make
- FastQC
- Picard
- GATK
- Samtools
- Python 2.7
- BWA
- Mummer (nucmer)
- LoFreq > 0.5
- pysam

Some scripts will complain if the dependencies cannot be found in
standard paths, but they should all allow you to point them to the
correct place with environment variables.


Usage
-----

ViPR is basically just the (Makefile based) glue for running a series
of scripts and external programs. The wrapper script is vipr.py

See

    $ vipr.py -h
    
for help.

Mandatory options are:
- -r `reference`: mapping reference fasta file
- -1 `reads`: first fastq[.gz] file (use -2 in addition for paired end sequencing)
- -p `PRIMER_FA`: a list of primers you used for amplification
- -o `OUT_DIR`: the output directory in which output files will be stored


Output
------

vipr.py produces a lot of files from coverage plots, to consensus
sequences and primer positions, as well as SNV predictions etc in the
specified output directory. The output directory will contain a
README.txt file that lists the important files. In brief, these are

- recal.bam: quality recalibrated mapping file with PCR duplicates removed
- mapping-success.txt: simple mapping success stats
- coverage.pdf: coverage plot
- cons_masked.fa: consensus/master sequence with primer regions masked
- final.snp: filtered, low-frequency SNV predictions



