#!/bin/bash


# Andreas Wilm <wilma@gis.a-star.edu.sg>
# Copyright: 2012, 2013 Genome Institute of Singapore
# License: GPL2


myname=$(basename $0)
echoinfo() {
    echo "INFO($myname): $@"
}
echowarn() {
    echo "WARNING($myname): $@" 1>&2
}
echoerror() {
    echo "ERROR($myname): $@" 1>&2
}
usage() {
    cat <<EOF
$myname: determine primer positions sequence from primer seqeunces

Options:
    -o | --output    main output file listing primer positions
    -p | --primers   primer fasta file
    -r | --ref       reference fasta file
       | --force     force overwriting of files
    -h | --help      print this help:
EOF
}

force=0
while [ "$1" != "" ]; do
    case $1 in
        -o | --output )
            shift
            output=$1
            ;;
        -p | --primers )
            shift
            primers_fa=$1
            ;;
        -r | --ref )
            shift
            refseq_fa=$1
            ;;
        --force )
            force=1
            ;;
        -h | --help ) 
            usage
            exit 0
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done

if [ ! $force -eq 1 ]; then
    # no accidential overwriting of files
    set -o noclobber
fi

if [ -z "$refseq_fa" ] && [ -z "$primers_fa" ] && [ -z "$output" ]; then
    echoerror "Required argument missing"
    echo
    usage
    exit 1
fi
for in_file in "$refseq_fa" "$primers_fa"; do
    if [ ! -s "$in_file" ]; then
        echoerror "input file $in_file does not exist or is empty. Exiting"
        exit 1
    fi
done


# setup output file names
#
primer_pos=${output}
delta=${output}.delta
coords=${output}.coords
nucmer_log=${output}.nucmer.log

# run nucmer
#
# reduce length and match criteria for nucmer
# nucmer won't run on windows/mac formatted primer files. try to convert and fail silently
#dos2unix $primers_fa >/dev/null 2>&1
# WARNING: may lead to race condition if this script is run concurrently!
file $primers_fa  | grep -q 'CRLF line terminator' && dos2unix $primers_fa >/dev/null 2>&1
if ! nucmer -c 10 -l 10 $refseq_fa $primers_fa -p $output > $nucmer_log 2>&1; then
    echoerror "nucmer failed. check $nucmer_log";
    exit 1
fi


# get coordinates
#
# -T tab-delimited
# -c percent coverage
# -d alignment direction
# -I min percent identity: -I 90 sometimes removes too much
# -l seq length info
# -r sort by ref-id and coords
show-coords  -T -c -d -r $delta  > $coords || exit 1


# convert to primer positions file
#
# note: not handling duplicates which shouldn't happen if primer design was not lousy
awk -v cov_thresh=80 '/^[0-9]/ {start=$1; end=$2; qcov=$9; ori=$11; qseq=$13;
  if (qcov>=cov_thresh) {printf "# %s\n", qseq; 
    if (ori==1) {printf "%d F\n", start} else {printf "%d R\n", end}}}' \
  $coords > $primer_pos || exit 1


# sanity check
#
num_primers=$(grep -c '^>' $primers_fa)
num_matches=$(grep -c '^[^#]' $primer_pos)
if [ $num_primers -ne $num_matches ]; then
    echowarn "found only $num_matches matches for a total of $num_primers primers" 1>&2
fi

echoinfo "main result file: $primer_pos"

exit 0

