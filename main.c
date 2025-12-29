/* main.c - Run runica_simple on EEG data */

#include "src/runica_simple.h"
#include "src/rt_nonfinite.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NCHANS 32
#define NPOINTS 30504
#define NDATA (NCHANS * NPOINTS)
#define NMATRIX (NCHANS * NCHANS)

int main(int argc, char *argv[])
{
    FILE *fp;
    float *data_float;
    double *data_double;
    const char *input_file;
    const char *weights_file;
    const char *sphere_file;
    int i;
    int nchan = NCHANS;
    int npoints = NPOINTS;
    boolean_T extended = true;  /* Default: extended ICA enabled */

    /* Parse command line arguments */
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input.fdt> <output.wts> <output.sph> [nchans] [npoints] [extended]\n", argv[0]);
        fprintf(stderr, "Example: %s data/eeglab_data.fdt data/eeglab_data.wts_darwin data/eeglab_data.sph_darwin 32 30504 1\n", argv[0]);
        fprintf(stderr, "\nOptional parameters:\n");
        fprintf(stderr, "  nchans:   Number of channels (default: %d)\n", NCHANS);
        fprintf(stderr, "  npoints:  Number of time points (default: %d)\n", NPOINTS);
        fprintf(stderr, "  extended: 1=extended ICA (tanh), 0=standard ICA (logistic) (default: 1)\n");
        return 1;
    }

    input_file = argv[1];
    weights_file = argv[2];
    sphere_file = argv[3];

    /* Optional nchan and npoints from command line */
    if (argc >= 5) {
        nchan = atoi(argv[4]);
    }
    if (argc >= 6) {
        npoints = atoi(argv[5]);
    }
    if (argc >= 7) {
        extended = (boolean_T)atoi(argv[6]);
    }

    /* Declare VLAs for weights and sphere matrices */
    const int nmatrix = nchan * nchan;
    const int ndata = nchan * npoints;
    double weights[nmatrix];
    double sphere[nmatrix];

    printf("======================================\n");
    printf("Running runica_simple\n");
    printf("======================================\n");
    printf("Input file: %s\n", input_file);
    printf("Channels: %d\n", nchan);
    printf("Data points: %d\n", npoints);
    printf("Mode: %s ICA\n", extended ? "Extended" : "Standard");
    printf("\n");

    /* Allocate memory */
    data_float = (float *)malloc(ndata * sizeof(float));
    data_double = (double *)malloc(ndata * sizeof(double));

    if (!data_float || !data_double) {
        fprintf(stderr, "Error: Failed to allocate memory\n");
        return 1;
    }

    /* Read input data (float32) */
    printf("Reading data file...\n");
    fp = fopen(input_file, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open input file: %s\n", input_file);
        free(data_float);
        free(data_double);
        return 1;
    }

    size_t nread = fread(data_float, sizeof(float), ndata, fp);
    fclose(fp);

    if (nread != (size_t)ndata) {
        fprintf(stderr, "Error: Expected %d values, read %zu\n", ndata, nread);
        free(data_float);
        free(data_double);
        return 1;
    }

    /* Convert float32 to float64 */
    printf("Converting to double precision...\n");
    for (i = 0; i < ndata; i++) {
        data_double[i] = (double)data_float[i];
    }
    free(data_float);

    /* Initialize and run ICA */
    printf("Initializing runica_simple...\n");
    runica_simple_initialize();

    printf("Running ICA...\n");
    runica_simple(data_double, weights, sphere, nchan, npoints, extended);

    printf("ICA completed successfully!\n\n");

    /* Save weights matrix */
    printf("Saving weights to: %s\n", weights_file);
    fp = fopen(weights_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write weights file: %s\n", weights_file);
        free(data_double);
        runica_simple_terminate();
        return 1;
    }
    fwrite(weights, sizeof(double), nmatrix, fp);
    fclose(fp);

    /* Save sphere matrix */
    printf("Saving sphere to: %s\n", sphere_file);
    fp = fopen(sphere_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write sphere file: %s\n", sphere_file);
        free(data_double);
        runica_simple_terminate();
        return 1;
    }

    /* Write sphere matrix (now real, not complex) */
    fwrite(sphere, sizeof(double), nmatrix, fp);
    fclose(fp);

    printf("\nDone!\n");
    printf("  Weights: %s (%dx%d matrix)\n", weights_file, nchan, nchan);
    printf("  Sphere:  %s (%dx%d matrix)\n", sphere_file, nchan, nchan);

    /* Cleanup */
    free(data_double);
    runica_simple_terminate();

    return 0;
}
