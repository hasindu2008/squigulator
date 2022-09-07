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

ex ./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -q a.fasta -n 10 --seed 1 --dwell-std 1.0 || die "Running the tool failed"
diff -q test/fasta.exp a.fasta || die "diff failed"
diff -q test/slow5.exp a.slow5 || die "diff failed"

ex ./sigsim -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o b.slow5 -q b.fastq -n 10 --seed 1 --prefix=yes || die "Running the tool failed"

echo "Test passed"