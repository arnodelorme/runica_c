# How to Create Parameterized runica_c

This document explains how to regenerate the C code with configurable parameters (extended, PCA, stop, maxsteps).

## Current Status

The existing `runica_simple.c` is generated from MATLAB with hardcoded parameters:
- Extended ICA: OFF (standard logistic)
- PCA: OFF
- maxsteps: 512
- stop: ~1e-6

## Steps to Create Parameterized Version

### 1. Create Parameterized MATLAB Function

In MATLAB, create a new file `runica_parameterized.m` based on `runica.m`:

```matlab
function [weights, sphere] = runica_parameterized(data, extended, pca, stop, maxsteps)
% RUNICA_PARAMETERIZED - Parameterized ICA for C code generation
%
% Inputs:
%   data      - [chans x frames] input data
%   extended  - 0=off, 1=on (extended ICA)
%   pca       - 0=off, N=keep N PCs
%   stop      - stopping threshold (e.g., 1e-6)
%   maxsteps  - maximum iterations (e.g., 512)
%
% Outputs:
%   weights   - [chans x chans] unmixing matrix
%   sphere    - [chans x chans] sphering matrix

% Set defaults
if nargin < 2, extended = 0; end
if nargin < 3, pca = 0; end
if nargin < 4, stop = 1e-6; end
if nargin < 5, maxsteps = 512; end

% Your ICA algorithm here using these parameters
% (Based on runica.m from EEGLAB)
...

end
```

### 2. Generate C Code with MATLAB Coder

```matlab
% In MATLAB with Coder toolbox installed:

% Define input types
data_type = coder.typeof(0, [32 Inf]);  % Variable time points
extended_type = coder.typeof(int32(0));
pca_type = coder.typeof(int32(0));
stop_type = coder.typeof(0.0);
maxsteps_type = coder.typeof(int32(0));

% Configure code generation
cfg = coder.config('lib');  % or 'exe' for executable
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

% Generate code
codegen -config cfg runica_parameterized ...
    -args {data_type, extended_type, pca_type, stop_type, maxsteps_type} ...
    -o runica_full
```

### 3. Update main.c to Accept Parameters

Modify `main.c` to parse command-line arguments:

```c
int main(int argc, char *argv[]) {
    // ...
    int extended = 0;
    int pca = 0;
    double stop = 1e-6;
    int maxsteps = 512;

    // Parse arguments
    if (argc > 4) extended = atoi(argv[4]);
    if (argc > 5) pca = atoi(argv[5]);
    if (argc > 6) stop = atof(argv[6]);
    if (argc > 7) maxsteps = atoi(argv[7]);

    // Call parameterized function
    runica_parameterized(data, extended, pca, stop, maxsteps,
                        weights, sphere);
    // ...
}
```

### 4. Create Test Matrix

Test all parameter combinations:

| Test | extended | pca | stop | maxsteps | Description |
|------|----------|-----|------|----------|-------------|
| 1    | 0        | 0   | 1e-6 | 512      | Standard ICA (default) |
| 2    | 1        | 0   | 1e-6 | 512      | Extended ICA |
| 3    | 0        | 16  | 1e-6 | 512      | PCA to 16 components |
| 4    | 1        | 16  | 1e-6 | 512      | Extended + PCA |
| 5    | 0        | 0   | 1e-7 | 512      | Tighter convergence |
| 6    | 0        | 0   | 1e-6 | 256      | Fewer iterations |

### 5. Validate Against MATLAB

For each test:
1. Run in C
2. Run in MATLAB with same parameters
3. Compare weights and sphere matrices
4. Max difference should be < 1e-10 (numerical precision)

## Alternative: Wrapper Approach

If full regeneration is not feasible, create separate compiled versions:

```bash
# Generate multiple versions
codegen ... -o runica_standard   # extended=0
codegen ... -o runica_extended   # extended=1
codegen ... -o runica_pca16      # pca=16
```

Then use a wrapper script to select the appropriate binary.

## Notes

- MATLAB Coder requires specific MATLAB versions and toolboxes
- Generated code may have platform-specific optimizations
- Always validate against MATLAB reference implementation
- Document exact MATLAB version and Coder version used
