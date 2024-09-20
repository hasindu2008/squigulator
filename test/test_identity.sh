#!/bin/bash

set -e

die() {
    echo "$@" >&2
    exit 1
}

REF_HG38=/genome/hg38noAlt.fa
REF_HG38_IDX=/genome/hg38noAlt.idx

REF_GENCODE=/genome/gencode.v40.transcripts.fa

test -e new.blow5 && rm new.blow5
test -e new.fastq && rm new.fastq

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

REMOVE_TMP(){
    rm -f new.blow5 new.fastq a.acc a.log
}

echo "R9 DNA ideal-time"
./squigulator -x dna-r9-prom --ideal-time $REF_HG38 -o new.blow5  2> a.log || die "squigulator failed"
eel  -i new.blow5 --config dna_r9.4.1_450bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identitydna.sh $REF_HG38_IDX new.fastq > a.acc  2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.97 a.acc
REMOVE_TMP

echo "R9 DNA"
./squigulator -x dna-r9-min $REF_HG38 -o new.blow5 2> a.log  || die "squigulator failed"
eel  -i new.blow5 --config dna_r9.4.1_450bps_sup.cfg --device cuda:all -o new.fastq &>> a.log || die "eel failed"
identitydna.sh $REF_HG38_IDX new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.95 a.acc
REMOVE_TMP

echo "R9 RNA"
./squigulator -x rna-r9-min $REF_GENCODE -o new.blow5 2> a.log || die "squigulator failed"
eel  -i new.blow5 --config rna_r9.4.1_70bps_hac.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF_GENCODE new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.75 a.acc
REMOVE_TMP

echo "R10 DNA ideal-time"
./squigulator -x dna-r10-min --ideal-time $REF_HG38 -o new.blow5 2> a.log || die "squigulator failed"
eel  -i new.blow5 --config dna_r10.4.1_e8.2_400bps_sup.cfg --device cuda:all -o new.fastq &>> a.log || die "eel failed"
identitydna.sh $REF_HG38_IDX new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.91 a.acc
REMOVE_TMP

echo "R10 DNA"
./squigulator -x dna-r10-prom $REF_HG38 -o new.blow5 -a new.sam 2> a.log || die "squigulator failed"
eel  -i new.blow5 --config dna_r10.4.1_e8.2_400bps_sup.cfg --device cuda:all -o new.fastq  &>> a.log || die "eel failed"
identitydna.sh $REF_HG38_IDX new.fastq > a.acc  2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.90 a.acc
samtools sort -o new.bam new.sam || die "samtools failed"
samtools index new.bam || die "samtools failed"
REMOVE_TMP

echo "RNA004 RNA"
./squigulator -x rna004-prom $REF_GENCODE -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF_GENCODE new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc
CHECK_ACC 0.77 a.acc
REMOVE_TMP

echo "CDNA R9"
./squigulator -x dna-r9-prom $REF_GENCODE -o new.blow5 --cdna 2> a.log || die "squigulator failed"
eel  -i new.blow5 --config dna_r9.4.1_450bps_sup.cfg --device cuda:all -o new.fastq &>> a.log || die "eel failed"
identitydna.sh $REF_GENCODE new.fastq > a.acc 2>> a.log || die "identitycdna failed"
cat a.acc
CHECK_ACC 0.94 a.acc
REMOVE_TMP



