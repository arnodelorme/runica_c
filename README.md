# runica_c

C implementation of runica for ICA decomposition of EEG data.

## Build

```bash
make
```

## Run

```bash
./run_runica.sh
```

Processes `data/eeglab_data.fdt` (32 channels, 30504 timepoints) and saves:
- `data/weights.bin` - unmixing matrix (32×32 float64)
- `data/sphere.bin` - sphering matrix (32×32 complex float64)

## Requirements

- macOS with Accelerate framework
- gcc compiler
