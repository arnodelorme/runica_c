#!/bin/bash
# run_all_tests.sh - Run all tests
# Usage: ./run_all_tests.sh [platform_suffix]
# Example: ./run_all_tests.sh _darwin

set -e

SUFFIX=${1:-"_test"}
TESTDIR="$(cd "$(dirname "$0")" && pwd)"

echo "========================================"
echo "Running All Tests"
echo "========================================"
echo "Platform suffix: $SUFFIX"
echo ""

PASSED=0
FAILED=0

# Test 1: Basic ICA
echo "----------------------------------------"
echo "Test 1: Basic ICA"
echo "----------------------------------------"
if bash "$TESTDIR/test_basic.sh" "$SUFFIX"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 2: Reproducibility
echo "----------------------------------------"
echo "Test 2: Reproducibility"
echo "----------------------------------------"
if bash "$TESTDIR/test_reproducibility.sh"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi
echo ""

# Summary
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All tests PASSED"
    exit 0
else
    echo "✗ Some tests FAILED"
    exit 1
fi
