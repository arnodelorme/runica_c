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
    creal_T sphere[NMATRIX];
    char weights_file[256];
    char sphere_file[256];
    int i;

    /* Default files */
    const char *input_file = "data/eeglab_data.fdt";
    const char *output_dir = "data";

    /* Allow command line override */
    if (argc > 1) {
        input_file = argv[1];
    }
    if (argc > 2) {
        output_dir = argv[2];
    }

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
    snprintf(weights_file, sizeof(weights_file), "%s/weights.bin", output_dir);
    printf("Saving weights to: %s\n", weights_file);
    fp = fopen(weights_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write weights file\n");
        free(data_double);
        runica_simple_terminate();
        return 1;
    }
    fwrite(weights, sizeof(double), NMATRIX, fp);
    fclose(fp);

    /* Save sphere matrix (real and imaginary parts) */
    snprintf(sphere_file, sizeof(sphere_file), "%s/sphere.bin", output_dir);
    printf("Saving sphere to: %s\n", sphere_file);
    fp = fopen(sphere_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write sphere file\n");
        free(data_double);
        runica_simple_terminate();
        return 1;
    }

    /* Write sphere as pairs of (real, imag) */
    for (i = 0; i < NMATRIX; i++) {
        fwrite(&sphere[i].re, sizeof(double), 1, fp);
        fwrite(&sphere[i].im, sizeof(double), 1, fp);
    }
    fclose(fp);

    printf("\nDone!\n");
    printf("  Weights: %s (%dx%d matrix)\n", weights_file, NCHANS, NCHANS);
    printf("  Sphere:  %s (%dx%d complex matrix)\n", sphere_file, NCHANS, NCHANS);

    /* Cleanup */
    free(data_double);
    runica_simple_terminate();

    return 0;
}
