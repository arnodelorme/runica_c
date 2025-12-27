#!/bin/bash
# run_plotica.sh - Plot ICA component topographies from runica output
# Usage: ./run_plotica.sh [dataset_basename] [n_channels] [suffix]
#
# Example: ./run_plotica.sh ./data/eeglab_data 32 _darwin

set -e  # Exit on error

# Default parameters
DATASET=${1:-"./data/eeglab_data"}
NCHANS=${2:-32}
SUFFIX=${3:-"_darwin"}

# Derived paths
SETFILE="${DATASET}.set"
DATADIR=$(dirname "$DATASET")
BASENAME=$(basename "$DATASET")
WTSFILE="${DATASET}.wts${SUFFIX}"
SPHFILE="${DATASET}.sph${SUFFIX}"

# MATLAB configuration
MATLAB_BIN="/Applications/MATLAB_R2025a.app/bin/matlab"
EEGLAB_PATH="~/eeglab"

echo "========================================"
echo "Plotting ICA topographies"
echo "========================================"
echo "Dataset: $DATASET"
echo "Channels: $NCHANS"
echo "Suffix: $SUFFIX"
echo ""

# Check if MATLAB is available
if [ ! -x "$MATLAB_BIN" ]; then
    echo "Error: MATLAB not found at $MATLAB_BIN"
    exit 1
fi

# Check if EEGLAB .set file exists
if [ ! -f "$SETFILE" ]; then
    echo "Error: EEGLAB .set file not found: $SETFILE"
    exit 1
fi

# Check if ICA files exist
if [ ! -f "$WTSFILE" ]; then
    echo "Error: Weights file not found: $WTSFILE"
    exit 1
fi

if [ ! -f "$SPHFILE" ]; then
    echo "Error: Sphere file not found: $SPHFILE"
    exit 1
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
pop_topoplot(EEG, 0, [1:${NCHANS}], 'ICA Component Topographies (runica_c)', [], 0, 'electrodes', 'off');

% Save outputs
pngfile = fullfile('${DATADIR_ABS}', '${BASENAME}_topoplot${SUFFIX}.png');
pdffile = fullfile('${DATADIR_ABS}', '${BASENAME}_topoplot${SUFFIX}.pdf');

print(gcf, pngfile, '-dpng', '-r150');
print(gcf, pdffile, '-dpdf', '-bestfit');

EEG = pop_saveset(EEG, 'filename', '${BASENAME}_with_ica${SUFFIX}.set', ...
    'filepath', '${DATADIR_ABS}/');

fprintf('Done! Saved:\\n');
fprintf('  %s\\n', pngfile);
fprintf('  %s\\n', pdffile);
fprintf('  ${DATADIR_ABS}/${BASENAME}_with_ica${SUFFIX}.set\\n');
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
    echo "Plotting completed!"
    echo "========================================"
    echo "Output files:"
    echo "  Topography PNG: ${DATASET}_topoplot${SUFFIX}.png"
    echo "  Topography PDF: ${DATASET}_topoplot${SUFFIX}.pdf"
    echo "  Dataset w/ICA:  ${DATADIR}/${BASENAME}_with_ica${SUFFIX}.set"
else
    echo "Error: MATLAB plotting failed"
    exit 1
fi

# Cleanup temporary files
rm -f "$MATLABSCRIPT"
