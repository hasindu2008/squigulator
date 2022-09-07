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

echo "Basic DNA"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -q a.fasta -n 10 --seed 1 --dwell-std 1.0 -r 20000 || die "Running the tool failed"
diff -q test/fasta.exp a.fasta || die "diff failed"
diff -q test/slow5.exp a.slow5 || die "diff failed"

echo "Basic RNA"
ex ./sigsim -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -q a.fastq -n 10 --seed 1 --prefix=yes --dwell-std 3.0 || die "Running the tool failed"
diff -q test/rna_slow5.exp a.slow5 || die "diff failed"

echo "--ideal"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal  -r 20000 || die "Running the tool failed"
diff -q test/dna_ideal_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-time"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-time  -r 20000 || die "Running the tool failed"
diff -q test/dna_ideal_time_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-amp"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-amp  -r 20000 || die "Running the tool failed"
diff -q test/dna_ideal_amp_slow5.exp a.slow5 || die "diff failed"

echo "--prefix=yes"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --prefix=yes  -r 20000 || die "Running the tool failed"
diff -q test/dna_prefix_slow5.exp a.slow5 || die "diff failed"

echo "--prefix=no"
ex ./sigsim -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -n 2 --seed 1 --dwell-std 3.0 || die "Running the tool failed"
diff -q test/rna_prefixno_slow5.exp a.slow5 || die "diff failed"

echo "--full-contigs"
ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 --seed 1 --full-contigs || die "Running the tool failed"
diff -q test/dna_full_contig.exp a.slow5 || die "diff failed"


echo "Test passed"