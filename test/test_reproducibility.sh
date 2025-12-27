#!/bin/bash
# test_reproducibility.sh - Verify deterministic results
# Runs ICA twice and compares outputs

set -e

TESTDIR="$(cd "$(dirname "$0")" && pwd)"
ROOTDIR="$(dirname "$TESTDIR")"

DATASET="../data/eeglab_data"
DATAFILE="${DATASET}.fdt"

WTS1="${TESTDIR}/repro_run1.wts"
SPH1="${TESTDIR}/repro_run1.sph"
WTS2="${TESTDIR}/repro_run2.wts"
SPH2="${TESTDIR}/repro_run2.sph"

echo "========================================"
echo "Reproducibility Test"
echo "========================================"
echo "Testing if ICA produces identical results on repeated runs"
echo ""

# Check if data file exists
if [ ! -f "$DATAFILE" ]; then
    echo "Error: Test data not found: $DATAFILE"
    exit 1
fi

# Find executable
RUNICA_EXE=""
for exe in runica_darwin runica_linux runica_expanse runica.exe runica; do
    if [ -f "$ROOTDIR/$exe" ]; then
        RUNICA_EXE="$ROOTDIR/$exe"
        break
    fi
done

if [ -z "$RUNICA_EXE" ]; then
    echo "Error: No runica executable found"
    exit 1
fi

# Run 1
echo "Run 1..."
$RUNICA_EXE "$DATAFILE" "$WTS1" "$SPH1" > /dev/null 2>&1

# Run 2
echo "Run 2..."
$RUNICA_EXE "$DATAFILE" "$WTS2" "$SPH2" > /dev/null 2>&1

# Compare using binary diff
echo ""
echo "Comparing outputs..."
if diff "$WTS1" "$WTS2" > /dev/null && diff "$SPH1" "$SPH2" > /dev/null; then
    echo "✓ Results are IDENTICAL (byte-for-byte)"
    echo "Test PASSED - ICA is deterministic"
    rm -f "$WTS1" "$SPH1" "$WTS2" "$SPH2"
    exit 0
else
    echo "✗ Results DIFFER"
    echo "Test FAILED - ICA is non-deterministic"
    echo "This is expected if random number generator is not seeded"
    exit 1
fi
