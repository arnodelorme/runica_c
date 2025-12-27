#!/bin/bash
# run_compica.sh - Compare ICA results from different platforms
# Usage: ./run_compica.sh [dataset_basename] [suffix1] [suffix2] [n_channels]
#
# Example: ./run_compica.sh ./data/eeglab_data _win _darwin 32

set -e  # Exit on error

# Default parameters
DATASET=${1:-"./data/eeglab_data"}
SUFFIX1=${2:-"_win"}
SUFFIX2=${3:-"_darwin"}
NCHANS=${4:-32}

# Derived paths
DATADIR=$(dirname "$DATASET")
BASENAME=$(basename "$DATASET")
WTS1="${DATASET}.wts${SUFFIX1}"
SPH1="${DATASET}.sph${SUFFIX1}"
WTS2="${DATASET}.wts${SUFFIX2}"
SPH2="${DATASET}.sph${SUFFIX2}"

# MATLAB configuration
MATLAB_BIN="/Applications/MATLAB_R2025a.app/bin/matlab"

echo "========================================"
echo "Comparing ICA results"
echo "========================================"
echo "Dataset: $DATASET"
echo "Platform 1: $SUFFIX1"
echo "Platform 2: $SUFFIX2"
echo "Channels: $NCHANS"
echo ""

# Check if MATLAB is available
if [ ! -x "$MATLAB_BIN" ]; then
    echo "Error: MATLAB not found at $MATLAB_BIN"
    exit 1
fi

# Check if files exist
for file in "$WTS1" "$SPH1" "$WTS2" "$SPH2"; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file"
        exit 1
    fi
done

# Create MATLAB script for comparison
MATLABSCRIPT="${DATASET}_compare.m"
DATADIR_ABS=$(cd "$DATADIR" && pwd)
echo "Creating MATLAB comparison script..."

cat > "$MATLABSCRIPT" << EOFMATLAB
% Auto-generated script to compare ICA results

% Load first set of matrices
fprintf('Loading ${SUFFIX1} matrices...\\n');
wts1file = fullfile('${DATADIR_ABS}', '${BASENAME}.wts${SUFFIX1}');
sph1file = fullfile('${DATADIR_ABS}', '${BASENAME}.sph${SUFFIX1}');

fid = fopen(wts1file, 'rb');
if fid == -1
    error('Cannot open weights file: %s', wts1file);
end
weights1 = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

fid = fopen(sph1file, 'rb');
if fid == -1
    error('Cannot open sphere file: %s', sph1file);
end
sphere1 = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

% Load second set of matrices
fprintf('Loading ${SUFFIX2} matrices...\\n');
wts2file = fullfile('${DATADIR_ABS}', '${BASENAME}.wts${SUFFIX2}');
sph2file = fullfile('${DATADIR_ABS}', '${BASENAME}.sph${SUFFIX2}');

fid = fopen(wts2file, 'rb');
if fid == -1
    error('Cannot open weights file: %s', wts2file);
end
weights2 = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

fid = fopen(sph2file, 'rb');
if fid == -1
    error('Cannot open sphere file: %s', sph2file);
end
sphere2 = fread(fid, [${NCHANS}, ${NCHANS}], 'float64')';
fclose(fid);

% Compare weights
fprintf('\\n======================================\\n');
fprintf('WEIGHTS COMPARISON\\n');
fprintf('======================================\\n');
wts_diff = weights1 - weights2;
wts_max_abs = max(abs(wts_diff(:)));
wts_max_rel = max(abs(wts_diff(:)) ./ (abs(weights1(:)) + eps));
fprintf('Max absolute difference: %.15e\\n', wts_max_abs);
fprintf('Max relative difference: %.15e\\n', wts_max_rel);
fprintf('Mean absolute difference: %.15e\\n', mean(abs(wts_diff(:))));
fprintf('Mean relative difference: %.15e\\n', mean(abs(wts_diff(:)) ./ (abs(weights1(:)) + eps)));

% Compare sphere
fprintf('\\n======================================\\n');
fprintf('SPHERE COMPARISON\\n');
fprintf('======================================\\n');
sph_diff = sphere1 - sphere2;
sph_max_abs = max(abs(sph_diff(:)));
sph_max_rel = max(abs(sph_diff(:)) ./ (abs(sphere1(:)) + eps));
fprintf('Max absolute difference: %.15e\\n', sph_max_abs);
fprintf('Max relative difference: %.15e\\n', sph_max_rel);
fprintf('Mean absolute difference: %.15e\\n', mean(abs(sph_diff(:))));
fprintf('Mean relative difference: %.15e\\n', mean(abs(sph_diff(:)) ./ (abs(sphere1(:)) + eps)));

% Overall summary
fprintf('\\n======================================\\n');
fprintf('SUMMARY\\n');
fprintf('======================================\\n');
if wts_max_abs < 1e-10 && sph_max_abs < 1e-10
    fprintf('Results are IDENTICAL (within machine precision)\\n');
elseif wts_max_rel < 1e-6 && sph_max_rel < 1e-6
    fprintf('Results are VERY CLOSE (relative diff < 1e-6)\\n');
elseif wts_max_rel < 1e-3 && sph_max_rel < 1e-3
    fprintf('Results are SIMILAR (relative diff < 1e-3)\\n');
else
    fprintf('Results DIFFER SIGNIFICANTLY\\n');
end
fprintf('\\n');
EOFMATLAB

echo "MATLAB script saved to: $MATLABSCRIPT"
echo ""

# Run MATLAB
echo "Running MATLAB comparison..."
FULLPATH=$(cd "$(dirname "${MATLABSCRIPT}")" && pwd)/$(basename "${MATLABSCRIPT}")
$MATLAB_BIN -batch "cd('$(pwd)'); run('${FULLPATH}')"

# Cleanup temporary files
rm -f "$MATLABSCRIPT"

echo ""
echo "Comparison complete!"
