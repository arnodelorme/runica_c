#!/bin/bash
# test_crossplatform.sh - Compare results across platforms
# Usage: ./test_crossplatform.sh [suffix1] [suffix2]
# Example: ./test_crossplatform.sh _darwin _win

set -e

TESTDIR="$(cd "$(dirname "$0")" && pwd)"
ROOTDIR="$(dirname "$TESTDIR")"

SUFFIX1=${1:-"_darwin"}
SUFFIX2=${2:-"_win"}

WTS1="${TESTDIR}/basic${SUFFIX1}.wts"
SPH1="${TESTDIR}/basic${SUFFIX1}.sph"
WTS2="${TESTDIR}/basic${SUFFIX2}.wts"
SPH2="${TESTDIR}/basic${SUFFIX2}.sph"

echo "========================================"
echo "Cross-Platform Comparison Test"
echo "========================================"
echo "Comparing: $SUFFIX1 vs $SUFFIX2"
echo ""

# Check if files exist
for file in "$WTS1" "$SPH1" "$WTS2" "$SPH2"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file not found: $file"
        echo ""
        echo "Please run test_basic.sh for both platforms first:"
        echo "  ./test_basic.sh $SUFFIX1"
        echo "  ./test_basic.sh $SUFFIX2"
        exit 1
    fi
done

# Create comparison using MATLAB
MATLAB_BIN="/Applications/MATLAB_R2025a.app/bin/matlab"

if [ ! -x "$MATLAB_BIN" ]; then
    echo "MATLAB not found. Comparing file sizes only..."
    ls -lh "$WTS1" "$WTS2" "$SPH1" "$SPH2"
    exit 0
fi

MATLABSCRIPT="${TESTDIR}/compare_temp.m"

cat > "$MATLABSCRIPT" << EOFMATLAB
% Compare test results

fprintf('Loading %s matrices...\\n', '${SUFFIX1}');
fid = fopen('${WTS1}', 'rb');
wts1 = fread(fid, [32, 32], 'float64')';
fclose(fid);
fid = fopen('${SPH1}', 'rb');
sph1 = fread(fid, [32, 32], 'float64')';
fclose(fid);

fprintf('Loading %s matrices...\\n', '${SUFFIX2}');
fid = fopen('${WTS2}', 'rb');
wts2 = fread(fid, [32, 32], 'float64')';
fclose(fid);
fid = fopen('${SPH2}', 'rb');
sph2 = fread(fid, [32, 32], 'float64')';
fclose(fid);

fprintf('\\nWeights comparison:\\n');
wts_diff = wts1 - wts2;
fprintf('  Max abs diff: %.3e\\n', max(abs(wts_diff(:))));
fprintf('  Max rel diff: %.3e\\n', max(abs(wts_diff(:)) ./ (abs(wts1(:)) + eps)));

fprintf('\\nSphere comparison:\\n');
sph_diff = sph1 - sph2;
fprintf('  Max abs diff: %.3e\\n', max(abs(sph_diff(:))));
fprintf('  Max rel diff: %.3e\\n', max(abs(sph_diff(:)) ./ (abs(sph1(:)) + eps)));

wts_ok = max(abs(wts_diff(:))) < 1e-10;
sph_ok = max(abs(sph_diff(:))) < 1e-10;

if wts_ok && sph_ok
    fprintf('\\n✓ Test PASSED - Results are identical\\n');
    exit(0);
else
    fprintf('\\n✗ Test FAILED - Results differ significantly\\n');
    exit(1);
end
EOFMATLAB

$MATLAB_BIN -batch "cd('${TESTDIR}'); run('compare_temp.m')"
RESULT=$?

rm -f "$MATLABSCRIPT"

if [ $RESULT -eq 0 ]; then
    echo "Cross-platform test PASSED"
else
    echo "Cross-platform test FAILED"
fi

exit $RESULT
