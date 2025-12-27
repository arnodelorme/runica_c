# runica_c

C implementation of runica for ICA decomposition of EEG data.

## Build

### Osx
```bash
make
```

### Windows
```bash
gcc -c src/runica_simple.c -o obj/runica_simple.o -O3 -Isrc
gcc -c src/rtGetInf.c -o obj/rtGetInf.o -O3 -Isrc
gcc -c src/rtGetNaN.c -o obj/rtGetNaN.o -O3 -Isrc
gcc -c src/rt_nonfinite.c -o obj/rt_nonfinite.o -O3 -Isrc
gcc -c main.c -o obj/main.o -O3 -Isrc
gcc -o runica.exe obj/runica_simple.o obj/rtGetInf.o obj/rtGetNaN.o obj/rt_nonfinite.o obj/main.o -lm
```

## Run

```bash
# Run ICA
./run_runica.sh [dataset] [nchans] [suffix]
./run_runica.sh ./data/eeglab_data 32 _darwin

# Plot topographies (requires MATLAB/EEGLAB)
./run_plotica.sh [dataset] [nchans] [suffix]

# Compare platforms
./run_compica.sh [dataset] [suffix1] [suffix2]
./run_compica.sh ./data/eeglab_data _darwin _win
```

Output files with platform suffix:
- `dataset.wts_suffix` - unmixing matrix (32×32 float64)
- `dataset.sph_suffix` - sphering matrix (32×32 float64)

## Current Implementation

**MATLAB Coder generated code** with hardcoded parameters:
- **Algorithm**: Standard Infomax ICA (logistic)
- **Extended ICA**: ❌ Not implemented (would require tanh for sub-Gaussian)
- **PCA**: ❌ Not implemented
- **maxsteps**: 512 (hardcoded)
- **stop**: ~1e-6 (hardcoded)
- **Channels**: 32 (hardcoded)
- **Random seed**: Fixed (deterministic results)

To enable configurable parameters (extended, PCA, stop, maxsteps), see `test/HOWTO_PARAMETERIZE.md`.

## Testing

```bash
cd test
./run_all_tests.sh _darwin    # Run all tests
./test_basic.sh _darwin        # Basic ICA test
./test_reproducibility.sh      # Verify deterministic results
./test_crossplatform.sh _darwin _win  # Compare platforms
```

See `test/README.md` for details.

## Requirements

- macOS with Accelerate framework (or Linux with BLAS)
- gcc compiler
- MATLAB with EEGLAB (for plotting and comparison)
