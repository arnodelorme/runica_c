#!/bin/bash
# run_runica.sh - Run runica_simple on EEG data and plot component topographies
# Usage: ./run_runica.sh [dataset_basename] [n_channels]
#
# Example: ./run_runica.sh ./data/eeglab_data 32

set -e  # Exit on error

# Detect platform
if [[ -n "$SLURM_CLUSTER_NAME" ]]; then
    PLATFORM="$SLURM_CLUSTER_NAME"
else
    PLATFORM="$(hostname -s)"
fi

echo "$PLATFORM"
if [ "$PLATFORM" = "MacBook-Pro-10" ]; then
    SUFFIX="_darwin"
else
    if [ "$PLATFORM" = "expanse" ]; then
        SUFFIX="_expanse"
    else
        SUFFIX="_linux"
    fi
fi

# Default parameters
DATASET=${1:-"./data/eeglab_data"}
NCHANS=${2:-32}

# Derived paths
DATAFILE="${DATASET}.fdt"
SETFILE="${DATASET}.set"
DATADIR=$(dirname "$DATASET")
BASENAME=$(basename "$DATASET")
WTSFILE="${DATASET}.wts${SUFFIX}"
SPHFILE="${DATASET}.sph${SUFFIX}"

# MATLAB configuration
MATLAB_BIN="/Applications/MATLAB_R2025a.app/bin/matlab"
EEGLAB_PATH="~/eeglab"

echo "========================================"
echo "Running runica_simple and plotting"
echo "========================================"
echo "Dataset: $DATASET"
echo "Channels: $NCHANS"
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
./runica_darwin "$DATAFILE" "$WTSFILE" "$SPHFILE"

if [ $? -ne 0 ]; then
    echo "Error: runica_darwin failed"
    exit 1
fi

echo ""
echo "ICA completed successfully!"
echo "  Weights: $WTSFILE"
echo "  Sphere:  $SPHFILE"
echo ""

# Only plot on Darwin
if [ "$PLATFORM" != "MacBook-Pro-10" ]; then
    echo "Skipping topography plotting (only available on macOS)"
    exit 0
fi

# Check if MATLAB and EEGLAB are available
if [ ! -x "$MATLAB_BIN" ]; then
    echo "Warning: MATLAB not found at $MATLAB_BIN"
    echo "Skipping topography plotting"
    exit 0
fi

if [ ! -f "$SETFILE" ]; then
    echo "Warning: EEGLAB .set file not found: $SETFILE"
    echo "Skipping topography plotting"
    echo "To create plots, ensure $SETFILE exists"
    exit 0
fi

# Create MATLAB script for plotting
MATLABSCRIPT="${DATASET}_plot.m"
DATADIR_ABS=$(cd "$DATADIR" && pwd)
echo "Creating MATLAB plotting script..."

cat > "$MATLABSCRIPT" << EOFMATLAB
% Auto-generated script to plot ICA topographies

addpath('$EEGLAB_PATH');
eeglab nogui;

% Load dataset
fprintf('Loading dataset...\\n');
EEG = pop_loadset('filename', '${BASENAME}.set', 'filepath', '${DATADIR_ABS}/');

% Load ICA matrices
fprintf('Loading ICA matrices...\\n');
wtsfile = fullfile('${DATADIR_ABS}', '${BASENAME}.wts${SUFFIX}');
sphfile = fullfile('${DATADIR_ABS}', '${BASENAME}.sph${SUFFIX}');

% Read weights matrix
fid = fopen(wtsfile, 'rb');
if fid == -1
    error('Cannot open weights file: %s', wtsfile);
end
weights = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

% Read sphere matrix (real)
fid = fopen(sphfile, 'rb');
if fid == -1
    error('Cannot open sphere file: %s', sphfile);
end
sphere = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

% Import into EEG structure
EEG.icaweights = weights;
EEG.icasphere = sphere;
EEG.icawinv = pinv(weights * sphere);
EEG.icaact = [];
EEG = eeg_checkset(EEG);

fprintf('Plotting component topographies...\\n');

% Plot all components
pop_topoplot(EEG, 0, [1:32], 'ICA Component Topographies (runica_c)', [], 0, 'electrodes', 'off');

% Save outputs
pngfile = fullfile('${DATADIR_ABS}', '${BASENAME}_topoplot.png');
pdffile = fullfile('${DATADIR_ABS}', '${BASENAME}_topoplot.pdf');

print(gcf, pngfile, '-dpng', '-r150');
print(gcf, pdffile, '-dpdf', '-bestfit');

EEG = pop_saveset(EEG, 'filename', '${BASENAME}_with_ica.set', ...
    'filepath', '${DATADIR_ABS}/');

fprintf('Done! Saved:\\n');
fprintf('  %s\\n', pngfile);
fprintf('  %s\\n', pdffile);
fprintf('  ${DATADIR_ABS}/${BASENAME}_with_ica.set\\n');
EOFMATLAB

echo "MATLAB script saved to: $MATLABSCRIPT"
echo ""

# Run MATLAB
echo "Running MATLAB to generate topography plots..."
FULLPATH=$(cd "$(dirname "${MATLABSCRIPT}")" && pwd)/$(basename "${MATLABSCRIPT}")
$MATLAB_BIN -batch "cd('$(pwd)'); run('${FULLPATH}')"

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "All done!"
    echo "========================================"
    echo "Output files:"
    echo "  ICA weights:    $WTSFILE"
    echo "  ICA sphere:     $SPHFILE"
    echo "  Topography PNG: ${DATASET}_topoplot.png"
    echo "  Topography PDF: ${DATASET}_topoplot.pdf"
    echo "  Dataset w/ICA:  ${DATADIR}/${BASENAME}_with_ica.set"
else
    echo "Warning: MATLAB plotting failed"
fi

# Cleanup temporary files
rm -f "$MATLABSCRIPT"
