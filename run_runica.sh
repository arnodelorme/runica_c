#!/bin/bash
# run_runica.sh - Run runica_simple on EEG data
# Usage: ./run_runica.sh [dataset_basename] [n_channels] [suffix]
#
# Example: ./run_runica.sh ./data/eeglab_data 32 _darwin

set -e  # Exit on error

# Default parameters
DATASET=${1:-"./data/eeglab_data"}
NCHANS=${2:-32}
SUFFIX=${3:-"_darwin"}
EXTENDED=${4:-1}  # 1=extended ICA (tanh), 0=standard ICA (logistic)

# Derived paths
DATAFILE="${DATASET}.fdt"
DATADIR=$(dirname "$DATASET")
BASENAME=$(basename "$DATASET")
WTSFILE="${DATASET}.wts${SUFFIX}"
SPHFILE="${DATASET}.sph${SUFFIX}"

echo "========================================"
echo "Running runica_simple"
echo "========================================"
echo "Dataset: $DATASET"
echo "Channels: $NCHANS"
echo "Suffix: $SUFFIX"
echo "Extended: $EXTENDED (1=extended, 0=standard)"
echo ""

# Check if data file exists
if [ ! -f "$DATAFILE" ]; then
    echo "Error: Data file $DATAFILE not found"
    exit 1
fi

# Find runica executable
RUNICA_EXE=""
for exe in runica_darwin runica_linux runica_expanse runica.exe runica; do
    if [ -f "./$exe" ]; then
        RUNICA_EXE="./$exe"
        break
    fi
done

if [ -z "$RUNICA_EXE" ]; then
    echo "No runica executable found. Building..."
    make clean
    make
    # Try again to find it
    for exe in runica_darwin runica_linux runica_expanse runica.exe runica; do
        if [ -f "./$exe" ]; then
            RUNICA_EXE="./$exe"
            break
        fi
    done
fi

if [ -z "$RUNICA_EXE" ]; then
    echo "Error: Could not find or build runica executable"
    exit 1
fi

# Run runica
# Note: executable uses defaults for npoints, but we pass nchans and extended
echo "Running $RUNICA_EXE..."
$RUNICA_EXE "$DATAFILE" "$WTSFILE" "$SPHFILE" "$NCHANS" 30504 "$EXTENDED"

if [ $? -ne 0 ]; then
    echo "Error: $RUNICA_EXE failed"
    exit 1
fi

echo ""
echo "========================================"
echo "ICA completed successfully!"
echo "========================================"
echo "Output files:"
echo "  Weights: $WTSFILE"
echo "  Sphere:  $SPHFILE"
echo ""
echo "To plot topographies, run:"
echo "  ./run_plotica.sh $DATASET $NCHANS $SUFFIX"
