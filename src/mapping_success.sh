#!/bin/bash


# Andreas Wilm <wilma@gis.a-star.edu.sg>
# Copyright: 2012, 2013 Genome Institute of Singapore
# License: GPL2

# Quick and dirty way to determine mapping success


usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): Determine mapping success (mapped reads in BAM vs reads in FastQ)

  Options:
    -b | --bam  <file>  : BAM file
    -f | --fastq <file> : Fastq file (if paired-end, simply one of them and use -p)
    -p | --pe           : This was a paired-end run
    -h | --help         : Display this help
EOF
}

pe=0
while [ "$1" != "" ]; do
    case $1 in
        -b | --bam )
            shift
            bam=$1
            ;;
        -f | --fastq )
            shift
            fastq=$1
            ;;
        -p | --pe )
            pe=1
            ;;
        -h | --help )
            usage
            exit
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done


# make sure all necessary args where given and files exist
#
if [ ! -e $fastq ] || [ -z $fastq ]; then
    echo "FATAL: fastq file \"$fastq\" missing" 1>&2
    echo
    usage
    exit 1
fi
if [ ! -e $bam ] || [ -z $bam ]; then
    echo "FATAL: BAM file \"$bam\" missing" 1>&2
    echo
    usage
    exit 1
fi


which samtools >/dev/null || exit 1

test -e ${bam}.bai || samtools index $bam

# if fastqc was run, used num reads from there
fastqcdata=${fastq%.gz}
fastqcdata=${fastqcdata%.fastq}
fastqcdata=${fastqcdata%.txt}
fastqcdata=${fastqcdata}_fastqc/fastqc_data.txt
if [ $fastqcdata -nt $fastq ]; then
    echo "Extracting number of reads from FastQC analysis ($fastqcdata) instead of FastQC (faster!)"
    nfastq=$(awk '/^Total S/ {print $NF; exit 1}' $fastqcdata)
else
    echo "Extracting number of reads from FastQ files (which might take some time...)"
    nfastq=$(zgrep -c '^@.*:' $fastq)
fi
if [ $pe -eq 1 ]; then
    nfastq=$(echo $nfastq*2 | bc)
fi

nmapped=$(samtools idxstats $bam | awk '{print $3; exit}')
    
echo  "$nmapped of $nfastq reads mapped"


