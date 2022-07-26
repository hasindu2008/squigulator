#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

./sigsim test/nCoV-2019.reference.fasta -o a.slow5 -q a.fasta -n 10 --seed 1 || die "Running the tool failed"
diff -q test/fasta.exp a.fasta || die "diff failed"
diff -q test/slow5.exp a.slow5 || die "diff failed"

echo "Test passed"