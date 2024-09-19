#!/bin/bash

set -e

die() {
    echo "$@" >&2
    exit 1
}

REF_HG38=/genome/hg38noAlt.fa
REF_HG38_IDX=/genome/hg38noAlt.idx

# should download the following from aws
METH_TRUTH=/home/hasindu/scratch/hg2_prom_lsk114_5khz/guppy_6.5.7_hac/PGXXXX230339_guppy657hac_mm226_f5c13_mfreq.tsv

MINIMOD=/data/hasindu/hasindu2008.git/minimod/minimod
COMPARE=~/hasindu2008.git/f5c/scripts/compare_methylation.py

REMOVE_TMP(){
    rm -f new.blow5 new.sam new.bam new.bam.bai new.bedmethyl methcomp.tsv a.acc a.log
}

REMOVE_TMP
test -e meth.tsv && rm meth.tsv
test -e ref_chr22.fa && rm ref_chr22.fa

tail -n +2 $METH_TRUTH | grep -w "chr22" | cut -f 1,2,7 > meth.tsv || die "failed extracting chr, pos, meth_freq"
samtools faidx ${REF_HG38} chr22 > ref_chr22.fa || die "failed extracting chr22 from ref"

CHECK_ACC(){
    THRESH=$1
    FILE=$2

    ACC=$(tail -1 $FILE | cut -f1)
    if [ $(echo "$ACC < $THRESH" | bc) -eq 1 ]; then
        die "FAILED: Accuracy $ACC is below threshold $THRESH"
    fi

    echo "PASSED: Accuracy $ACC is above threshold $THRESH"
    echo "________________________________________________"
    echo ""
    echo ""
}



RUN_TEST(){
    PROF=$1
    MODEL=$2
    ./squigulator ref_chr22.fa -x ${PROF} -f 10 -t 20 -o new.blow5 --meth-freq meth.tsv 2> a.log || die "squigulator failed"
    eel  -i new.blow5 --config ${MODEL} --device cuda:all -o new.sam --call_mods &>> a.log || die "eel failed"
    samtools fastq -TMM,ML new.sam  2>> a.log  | minimap2 -ax map-ont -y -Y --secondary=no ref_chr22.fa -  2>> a.log | samtools sort - -o new.bam  2>> a.log || die "samtools failed"
    samtools index new.bam || die "samtools index failed"
    # /install/modkit-v0.1.13/modkit pileup --cpg --ref ref_chr22.fa --ignore h -t 32 new.bam  new.tmp.bedmethyl
    # grep "chr22" new.tmp.bedmethyl |  grep -v nan  > new.bedmethyl
    ${MINIMOD} mod-freq ref_chr22.fa new.bam -b > new.bedmethyl 2>> a.log || die "minimod failed"
    ${COMPARE} ${METH_TRUTH} new.bedmethyl > methcomp.tsv 2>> a.log || die "compare failed"
    tail -n+2 methcomp.tsv | datamash ppearson 3:5 > a.acc || die "pearson failed"
    # ~/hasindu2008.git/f5c/scripts/plot_methylation.R  -i methcomp.tsv -o methcomp.pdf
    cat a.acc
    CHECK_ACC 0.93 a.acc
    REMOVE_TMP
}

echo "R9 DNA methylation"
RUN_TEST "dna-r9-prom" "dna_r9.4.1_450bps_modbases_5mc_cg_hac_prom.cfg"

# echo "R10 DNA methylation"
# RUN_TEST "dna-r10-prom" "dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_prom.cfg"

