# runica_c Test Suite

This directory contains tests for the runica_c implementation.

## Current Implementation Status

The current `runica_simple.c` is MATLAB Coder generated code with **hardcoded parameters**:

- **Extended ICA**: ❌ Not implemented (standard logistic ICA only)
- **PCA**: ❌ Not implemented
- **maxsteps**: ✅ Hardcoded to 512
- **stop**: ✅ Hardcoded convergence threshold (~1e-6)
- **lrate**: ✅ Hardcoded learning rate (~5e-4)

## To Enable Configurable Parameters

To add parameters (extended, PCA, stop, maxsteps), you need to:

1. Modify the MATLAB source: `~/eeglab/functions/sigprocfunc/runica.m`
2. Create a parameterized version that accepts these as inputs
3. Regenerate C code using MATLAB Coder with these parameters as codegen inputs
4. Replace the generated `runica_simple.c` and `runica_simple.h`

Example MATLAB Coder command structure:
```matlab
% In MATLAB with Coder:
cfg = coder.config('exe');
cfg.EnableOpenMP = false;
codegen -config cfg runica_parameterized ...
    -args {coder.typeof(0,[32 Inf]), ... % data
           coder.Constant(1), ...         % extended (0 or 1)
           coder.Constant(0), ...         % pca (0 = off)
           coder.Constant(1e-6), ...      % stop
           coder.Constant(512)}           % maxsteps
```

## Current Tests

The tests below verify the current implementation and provide a framework for cross-platform comparison:

- `test_basic.sh` - Basic ICA on test data
- `test_reproducibility.sh` - Verify deterministic results (same random seed)
- `test_crossplatform.sh` - Compare results across platforms

## Running Tests

```bash
cd test
./run_all_tests.sh
```

## Cross-Platform Comparison

Save results with platform suffix for comparison:
```bash
# On each platform:
./test_basic.sh _darwin    # macOS
./test_basic.sh _linux     # Linux
./test_basic.sh _win       # Windows

# Compare:
../run_compica.sh ./test_data _darwin _linux
```
