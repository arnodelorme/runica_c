#!/bin/bash
# run_runica.sh - Run runica_simple on EEG data
# Usage: ./run_runica.sh [dataset_basename]
#
# Example: ./run_runica.sh ./data/eeglab_data

set -e  # Exit on error

# Default parameters
DATASET=${1:-"./data/eeglab_data"}

# Derived paths
DATAFILE="${DATASET}.fdt"
DATADIR=$(dirname "$DATASET")

echo "========================================"
echo "Running runica_simple"
echo "========================================"
echo "Dataset: $DATASET"
echo ""

# Check if data file exists
if [ ! -f "$DATAFILE" ]; then
    echo "Error: Data file $DATAFILE not found"
    exit 1
fi

# Check if executable exists
if [ ! -f "./runica_darwin" ]; then
    echo "Executable not found. Building..."
    make clean
    make
    echo ""
fi

# Run runica
echo "Running runica_darwin..."
./runica_darwin "$DATAFILE" "$DATADIR"

if [ $? -ne 0 ]; then
    echo "Error: runica_darwin failed"
    exit 1
fi

echo ""
echo "========================================"
echo "All done!"
echo "========================================"
echo "Output files:"
echo "  Weights: ${DATADIR}/weights.bin"
echo "  Sphere:  ${DATADIR}/sphere.bin"
