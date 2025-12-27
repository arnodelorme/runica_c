#!/bin/bash
# test_basic.sh - Basic ICA test on eeglab data
# Usage: ./test_basic.sh [suffix]
# Example: ./test_basic.sh _darwin

set -e

SUFFIX=${1:-"_test"}
TESTDIR="$(cd "$(dirname "$0")" && pwd)"
ROOTDIR="$(dirname "$TESTDIR")"

# Use the main data file
DATASET="../data/eeglab_data"
DATAFILE="${DATASET}.fdt"
WTSFILE="${TESTDIR}/basic${SUFFIX}.wts"
SPHFILE="${TESTDIR}/basic${SUFFIX}.sph"

echo "========================================"
echo "Basic ICA Test"
echo "========================================"
echo "Input: $DATAFILE"
echo "Output suffix: $SUFFIX"
echo ""

# Check if data file exists
if [ ! -f "$DATAFILE" ]; then
    echo "Error: Test data not found: $DATAFILE"
    echo "Please ensure data/eeglab_data.fdt exists"
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

# Run ICA
echo "Running: $RUNICA_EXE"
$RUNICA_EXE "$DATAFILE" "$WTSFILE" "$SPHFILE"

if [ $? -eq 0 ]; then
    echo ""
    echo "Test PASSED"
    echo "Output files:"
    echo "  Weights: $WTSFILE"
    echo "  Sphere:  $SPHFILE"
    ls -lh "$WTSFILE" "$SPHFILE"
else
    echo "Test FAILED"
    exit 1
fi
