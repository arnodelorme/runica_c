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
./run_runica.sh
```

Processes `data/eeglab_data.fdt` (32 channels, 30504 timepoints) and saves:
- `data/weights.bin` - unmixing matrix (32×32 float64)
- `data/sphere.bin` - sphering matrix (32×32 complex float64)

## Requirements

- macOS with Accelerate framework
- gcc compiler
