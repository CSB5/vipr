ViPR
====

This is a collection of scripts used for the analysis of RNA viruses.
The main wrapper is vipr.py ("viper") which will create a makefile
defining the whole pipeline and subsequently call make.

Some of these scripts are quite specific to small genomes, others you
will able to reuse some of them for larger - e.g. bacterial - genomes.

WARNING: this pipeline doesn't use a framework, has a huge dependency list
and might make assumptions met only in our computing environment. Use at
your own risk.

Installation
------------

The impatient can simply use the full path for calling the scripts in
./src or add that directory to their PATH variable.

For a proper installation, cd into this directory and then use
$ python setup.py install --prefix
for installing the scripts.

For more option see
$ python setup.py install -h

Note, if you're using the prefix argument make sure that the
corresponding bin directory is in PATH and the corresponding python
library path is in PYTHONPATH.


Dependencies
------------

A (very likely) incomplete list of dependencies to be installed:
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


