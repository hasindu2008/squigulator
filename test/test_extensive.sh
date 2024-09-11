#!/bin/bash

set -e

die() {
    echo "$@" >&2
    exit 1
}

test/test_identity.sh || die "test_identity.sh failed"

echho "all done"