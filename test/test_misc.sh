#!/bin/bash

set -e

die() {
    echo "$@" >&2
    exit 1
}


REF_GENCODE=/genome/gencode.v40.transcripts.fa

REF_TRANS_BAM=/home/hasindu/scratch/uhr_prom_rna004/PNXRXX240011_reads_500k.bam

REMOVE_TMP(){
    rm -f new.blow5 new.fastq a.acc a.log a.paf count_joined.tsv count_new.tsv
}

REMOVE_TMP
test -e count.tsv && rm count.tsv


CHECK_ACC(){
    THRESH=$1
    FILE=$2

    ACC=$(tail -1 $FILE | cut -f2)
    if [ $(echo "$ACC < $THRESH" | bc) -eq 1 ]; then
        die "FAILED: Accuracy $ACC is below threshold $THRESH"
    fi

    echo "PASSED: Accuracy $ACC is above threshold $THRESH"
    echo "________________________________________________"
    echo ""
    echo ""
}



samtools  view -F 4 ${REF_TRANS_BAM} | cut -f 3 | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1 > count.tsv || die "failed extracting chr, pos, meth_freq"

RUN_REST(){
    minimap2 -cx map-ont --secondary=no  ${REF_GENCODE} new.fastq > a.paf  2>> a.log
    cat a.paf | cut -f 6 | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' > count_new.tsv || die "failed getting counts"
    join count.tsv count_new.tsv > count_joined.tsv || die "failed joining counts"
    cat count_joined.tsv | awk '{print $1"\t"$2"\t"$3}' | datamash ppearson 2:3 > a.acc || die "pearson failed"
    cat a.acc
    CHECK_ACC 0.95 a.acc
    REMOVE_TMP
}

PROF=rna004-prom
MODEL=rna_rp4_130bps_sup.cfg
./squigulator ${REF_GENCODE} -x ${PROF}  -t 20 -n 50000 -o new.blow5 --trans-count count.tsv 2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config ${MODEL} --device cuda:all -o new.fastq  &>> a.log || die "eel failed"
RUN_REST

PROF=dna-r9-prom
MODEL=dna_r9.4.1_450bps_hac_prom.cfg
./squigulator ${REF_GENCODE} -x ${PROF} -t 20 -n 50000 -o new.blow5 --trans-count count.tsv --cdna 2> a.log || die "squigulator failed"
eel  -i new.blow5 --config ${MODEL} --device cuda:all -o new.fastq  &>> a.log || die "eel failed"
RUN_REST