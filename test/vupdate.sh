#!/bin/bash

die() {
    echo "$@" >&2
    exit 1
}

[ $# -eq 2 ] || die "Usage: $0 <old_version> <new_version>"

OLD=$1
NEW=$2

grep "VN:$OLD" test/rna_paf.sam.exp || die "Version $OLD not found in test/rna_paf.sam.exp"
grep "VN:$OLD" test/dna_r10_paf-ref.sam.exp || die "Version $OLD not found in dna_r10_paf-ref.sam.exp"
grep "SQ_VERSION \"$OLD\"" src/version.h || die "Version $OLD not found in src/sq.h"
grep "VERSION=$OLD" README.md || die "Version $OLD not found in README.md"

sed -i "s/VN:$OLD/VN:$NEW/" test/rna_paf.sam.exp | grep "VN:$NEW"
sed -i "s/VN:$OLD/VN:$NEW/" test/dna_r10_paf-ref.sam.exp | grep "VN:$NEW"
sed -i "s/SQ_VERSION \"$OLD\"/SQ_VERSION \"$NEW\"/" src/version.h | grep "SQ_VERSION \"$NEW\""
sed -i "s/VERSION=$OLD/VERSION=$NEW/" README.md | grep "VERSION=$NEW"

grep "VN:$NEW" test/rna_paf.sam.exp || die "updated version $NEW not found in test/rna_paf.sam.exp"
grep "VN:$NEW" test/dna_r10_paf-ref.sam.exp || die "updated Version $NEW not found in dna_r10_paf-ref.sam.exp"
grep "SQ_VERSION \"$NEW\"" src/version.h || die "updated Version $NEW not found in src/sq.h"
grep "VERSION=$NEW" README.md || die "updated Version $NEW not found in README.md"


