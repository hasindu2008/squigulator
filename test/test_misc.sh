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




