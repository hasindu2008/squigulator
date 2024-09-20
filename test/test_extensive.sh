#!/bin/bash

set -e

die() {
    echo "$@" >&2
    exit 1
}

test/test_identity.sh || die "test_identity.sh failed"
test/test_meth.sh || die "test_meth.sh failed"
test/test_misc.sh || die "test_misc.sh failed"

echho "all done"