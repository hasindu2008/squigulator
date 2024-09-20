#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

# DNA R9
echo "Basic DNA"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -q a.fasta -n 10 --seed 1 --dwell-std 1.0 -r 20000 -t1 || die "Running the tool failed"
diff -q test/fasta.exp a.fasta || die "diff failed"
diff -q test/slow5.exp a.slow5 || die "diff failed"

# RNA R9
echo "Basic RNA"
ex ./squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -q a.fastq -n 10 --seed 1 --prefix=yes --dwell-std 3.0 -t1  || die "Running the tool failed"
diff -q test/rna_slow5.exp a.slow5 || die "diff failed"

# --ideal
echo "--ideal"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal  -r 20000 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-time"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-time  -r 20000 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_time_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-amp"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-amp  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_amp_slow5.exp a.slow5 || die "diff failed"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --amp-noise 0.0  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_amp_slow5.exp a.slow5 || die "diff failed"

# prefixes for dna
echo "--prefix=yes"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --prefix=yes  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_prefix_slow5.exp a.slow5 || die "diff failed"

# prefixes for rna
echo "--prefix=yes"
ex ./squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -n 2 --seed 1 --dwell-std 3.0 -t1 --prefix=yes || die "Running the tool failed"
diff -q test/rna_prefixyes_slow5.exp a.slow5 || die "diff failed"

echo "--prefix=no"
ex ./squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -n 2 --seed 1 --dwell-std 3.0 -t1 || die "Running the tool failed"
diff -q test/rna_prefixno_slow5.exp a.slow5 || die "diff failed"

# full contigs
echo "--full-contigs"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 --seed 1 --full-contigs  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_full_contig.exp a.slow5 || die "diff failed"

# r10 and paf and fasta
echo "r10 PAF out"
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 1 --seed 1 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf -q a.fa
diff -q test/dna_r10_paf.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf.paf.exp a.paf || die "diff failed"
diff -q test/dna_r10_paf.fa.exp a.fa || die "diff failed"

# r9 rna paf and sam and fasta outputs
echo "r9 rna paf out and sam out"
ex ./squigulator -x rna-r9-prom -o a.slow5 -n 1 --seed 1 --dwell-std 3.0 -t1 -t1 test/rnasequin_sequences_2.4.fa -c a.paf -q a.fa -a a.sam
diff -q test/rna_paf.exp a.slow5 || die "diff failed"
diff -q test/rna_paf.paf.exp a.paf || die "diff failed"
diff -q test/rna_paf.sam.exp a.sam || die "diff failed"
diff -q test/rna_paf.fa.exp a.fa || die "diff failed"

# --paf-ref and samout
echo "r10 dna paf out --paf-ref"
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 2 --seed 2 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf --paf-ref -a a.sam
diff -q test/dna_r10_paf-ref.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf-ref.paf.exp a.paf || die "diff failed"
diff -q test/dna_r10_paf-ref.sam.exp a.sam || die "diff failed"

# samout only
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 2 --seed 2 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf -a a.sam
diff -q test/dna_r10_paf-ref.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf-ref.sam.exp a.sam || die "diff failed"

# rna004
echo "rna004 test"
ex ./squigulator -x rna004-prom -o a.slow5 -n 1 --seed 1 --dwell-std 3.0 -t1 test/rnasequin_sequences_2.4.fa
diff -q test/rna004.slow5.exp a.slow5 || die "diff failed"

# -r and -f and --amp-noise
echo "read lens and fold coverage"
ex ./squigulator -x dna-r10-prom -o a.slow5 -r 20000 -f 1 --seed 2 --amp-noise 0.5 -t1 test/nCoV-2019.reference.fasta
diff -q test/dna_r10_amp_noise.exp a.slow5 || die "diff failed"

# dwell mean and std
ex ./squigulator -x rna004-min -o a.slow5 -n 1 --seed 1 --dwell-mean 30 --dwell-std 3.0 -t1 test/rnasequin_sequences_2.4.fa
diff -q test/rna004_dwell.exp a.slow5 || die "diff failed"

# bps
ex ./squigulator -x dna-r10-prom -o a.slow5 --seed 1 --bps 200 -t1 -n 2 test/nCoV-2019.reference.fasta
diff -q test/bps.exp a.slow5  || die "diff failed"

# cdna
ex ./squigulator -x dna-r10-min -o a.slow5 -n 1 --seed 1 --dwell-std 3.0 -t1 test/rnasequin_sequences_2.4.fa --cdna
diff -q test/cdna.exp a.slow5 || die "diff failed"

# trans count
ex ./squigulator -x rna004-prom -o a.slow5 -n 3 --seed 3 --trans-count test/sequin_count.tsv -t1 test/rnasequin_sequences_2.4.fa
diff -q test/trans_count.exp a.slow5 || die "diff failed"

# trans count cdna
ex ./squigulator -x dna-r10-min -o a.slow5 -n 3 --seed 3 --trans-count test/sequin_count.tsv -t1 test/rnasequin_sequences_2.4.fa --cdna
diff -q test/trans_count_cdna.exp a.slow5 || die "diff failed"

# trans trunc
ex ./squigulator -x rna004-prom -o a.slow5 -n 1 --seed 1 --trans-trunc -t1 test/rnasequin_sequences_2.4.fa
diff -q test/trans_trunc.exp a.slow5 || die "diff failed"

# ont-friendly
ex ./squigulator -x dna-r10-min -o a.slow5 -n 1 --seed 1 -t1 test/rnasequin_sequences_2.4.fa --ont-friendly=yes
diff -q test/ont_friendly.exp a.slow5 || die "diff failed"

# digitisation, sample-rate, range, offset-mean, offset-std, median-before-mean, median-before-std
ex ./squigulator -x dna-r10-min -o a.slow5 -n 1 --seed 1 -t1 test/rnasequin_sequences_2.4.fa --digitisation 4096 --sample-rate 10000 --range 300 --offset-mean -1000 --offset-std 0 --median-before-mean 100 --median-before-std 0
diff -q test/dev.exp a.slow5 || die "diff failed"

#meth r9
ex ./squigulator -x dna-r9-prom -o a.slow5 --seed 1 -t1 -n 2 -r 29000 test/nCoV-2019.reference.fasta --meth-freq test/mfreq.tsv
diff -q test/r9_mfreq.exp a.slow5  || die "diff failed"
# /install/buttery-eel-0.3.1+6.5.7/scripts/eel -i a.slow5 -o a.fastq --config  dna_r9.4.1_450bps_sup.cfg -x cuda:all
# minimap2 -ax map-ont /genome/nCoV-2019.reference.fasta  a.fastq --secondary=no | samtools sort - -o a.bam && samtools index a.bam
# f5c index a.fastq --slow5 a.slow5 && f5c call-methylation -r a.fastq  -g test/nCoV-2019.reference.fasta -b a.bam --slow5 a.slow5 -o meth.tsv && f5c meth-freq -i meth.tsv -s -o meth-freq.tsv

#meth r10
# ex ./squigulator -x dna-r10-prom -o a.slow5 --seed 1 -t1 -n 2 -r 29000 test/nCoV-2019.reference.fasta --meth-freq test/methfreq.tsv
# diff -q test/r10_methfreq.exp a.slow5  || die "diff failed"
# /install/buttery-eel-0.3.1+6.5.7/scripts/eel -i a.slow5 -o a.fastq --config  dna_r10.4.1_e8.2_400bps_5khz_sup.cfg -x cuda:all
# minimap2 -ax map-ont /genome/nCoV-2019.reference.fasta  a.fastq --secondary=no | samtools sort - -o a.bam && samtools index a.bam
# f5c index a.fastq --slow5 a.slow5 && f5c call-methylation -r a.fastq  -g test/nCoV-2019.reference.fasta -b a.bam --slow5 a.slow5 -o meth.tsv && f5c meth-freq -i meth.tsv -s -o meth-freq.tsv

# threads and batch size
redundancy_check () {
    N=$(grep -v ^[@#] a.slow5 | cut -f ${1}  | sort | uniq -c | sort -nr -k1,1 | head -1 | awk '{print $1}')
    [ "$N" != "1" ] && die "failed thread test for column ${1}"
}

ex ./squigulator -x dna-r9-min test/nCoV-2019.reference.fasta -n 100 -t 8 -K 10 -o a.slow5 --seed 1
# read_id        read_group      digitisation    offset  range   sampling_rate   len_raw_signal  raw_signal
# 9channel_number  10median_before   11read_number     12start_mux        13start_time
redundancy_check 1
redundancy_check 4
redundancy_check 7
redundancy_check 8
redundancy_check 10
redundancy_check 11
redundancy_check 13

echo "Test passed"