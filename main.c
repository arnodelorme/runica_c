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
    double weights[NMATRIX];
    double sphere[NMATRIX];
    const char *input_file;
    const char *weights_file;
    const char *sphere_file;
    int i;

    /* Parse command line arguments */
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input.fdt> <output.wts> <output.sph>\n", argv[0]);
        fprintf(stderr, "Example: %s data/eeglab_data.fdt data/eeglab_data.wts_darwin data/eeglab_data.sph_darwin\n", argv[0]);
        return 1;
    }

    input_file = argv[1];
    weights_file = argv[2];
    sphere_file = argv[3];

    printf("======================================\n");
    printf("Running runica_simple\n");
    printf("======================================\n");
    printf("Input file: %s\n", input_file);
    printf("Channels: %d\n", NCHANS);
    printf("Data points: %d\n", NPOINTS);
    printf("\n");

    /* Allocate memory */
    data_float = (float *)malloc(NDATA * sizeof(float));
    data_double = (double *)malloc(NDATA * sizeof(double));

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

    size_t nread = fread(data_float, sizeof(float), NDATA, fp);
    fclose(fp);

    if (nread != NDATA) {
        fprintf(stderr, "Error: Expected %d values, read %zu\n", NDATA, nread);
        free(data_float);
        free(data_double);
        return 1;
    }

    /* Convert float32 to float64 */
    printf("Converting to double precision...\n");
    for (i = 0; i < NDATA; i++) {
        data_double[i] = (double)data_float[i];
    }
    free(data_float);

    /* Initialize and run ICA */
    printf("Initializing runica_simple...\n");
    runica_simple_initialize();

    printf("Running ICA...\n");
    runica_simple(data_double, weights, sphere);

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
    fwrite(weights, sizeof(double), NMATRIX, fp);
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
    fwrite(sphere, sizeof(double), NMATRIX, fp);
    fclose(fp);

    printf("\nDone!\n");
    printf("  Weights: %s (%dx%d matrix)\n", weights_file, NCHANS, NCHANS);
    printf("  Sphere:  %s (%dx%d matrix)\n", sphere_file, NCHANS, NCHANS);

    /* Cleanup */
    free(data_double);
    runica_simple_terminate();

    return 0;
}
