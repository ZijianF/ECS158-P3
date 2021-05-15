#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>
double limit;
double complex p;

unsigned char *mandelbrot_mapping(size_t N, size_t C, double inc, double x_min, double y_max)
{
    p = CMPLX(2, 0);
    limit = 2.0;
    unsigned char *g_map = aligned_alloc(64, N * N * sizeof(unsigned char));
    unsigned char iteration = (unsigned char)C;
    unsigned char cutoff = 0;
    memset(g_map, 0, N * N * sizeof(unsigned char));

    for (size_t y = 0; y < N; y++)
    {
        double cpl_c = y_max - y * inc;
        for (size_t x = 0; x < N; x++)
        {
            double real_c = x_min + x * inc;
            double complex c = CMPLX(real_c, cpl_c);
            double complex zn = CMPLX(0, 0);
            cutoff = iteration;
            while (cutoff > 0)
            {

                double zn_abs = sqrt(creal(zn) * creal(zn) + cimag(zn) * cimag(zn));
                if (zn_abs >= limit)
                    break;
                zn = zn * zn + c;
                cutoff--;
            }
            g_map[y * N + x] = iteration - cutoff;
        }
    }

    return g_map;
}

// Write pgm format
void write_image(unsigned char *g_map, size_t N, double x_center, double y_center, double zoom, size_t iteration)
{
    char *fn = calloc(100, sizeof(char));
    char *info = calloc(100, sizeof(char));
    sprintf(fn, "mandel_%zu_%.3lf_%.3lf_%.3lf_%zu.pgm", N, x_center, y_center, zoom, iteration);
    sprintf(info, "P5\n%lu %lu\n%lu\n", N, N, iteration);
    FILE *pgm_fd = fopen(fn, "w");
    fwrite(info, sizeof(char), strlen(info), pgm_fd);

    fwrite(g_map, sizeof(unsigned char), N * N, pgm_fd);
    fclose(pgm_fd);
    free(g_map);
    free(fn);
}

int main(int argc, char *argv[])
{
    int N, C; // order, cutoff
    double x_center, y_center, zoom;
    if (argc < 6)
    {
        fprintf(stderr,
                "Usage: ./mandelbrot_serial order xcenter ycenter zoom cutoff\n");
        exit(1);
    }
    sscanf(argv[1], " %d", &N);
    sscanf(argv[2], " %lf ", &x_center);
    sscanf(argv[3], " %lf ", &y_center);
    sscanf(argv[4], " %lf ", &zoom);
    sscanf(argv[5], " %d ", &C);
    if (N < 128 || N > 8192)
    {
        fprintf(stderr, "Error: wrong order (128 <= N <= 8192)\n");
        exit(1);
    }
    if (x_center < -10.0 || x_center > 10.0)
    {
        fprintf(stderr, "Error: wrong x-center (-10.0 <= N <= 10.0)\n");
        exit(1);
    }
    if (y_center < -10.0 || y_center > 10.0)
    {
        fprintf(stderr, "Error: wrong y-center (-10.0 <= N <= 10.0)\n");
        exit(1);
    }
    if (zoom < 0.0 || zoom > 100.0)
    {
        fprintf(stderr, "Error: wrong zoom (0.0 <= N <= 100.0)\n");
        exit(1);
    }
    if (C < 10 || C > 255)
    {
        fprintf(stderr, "Error: wrong cutoff (10 <= N <= 255)\n");
        exit(1);
    }
    double complex *map = aligned_alloc(64, N * N * sizeof(double complex));
    memset(map, 0, N * N * sizeof(double complex));

    double inc = pow(2, (-1) * zoom);
    double width = inc * N;
    double x_min = x_center - width / 2;
    double y_max = y_center + width / 2;

    unsigned char *gray_map = mandelbrot_mapping(N, C, inc, x_min, y_max);

    write_image(gray_map, N, x_center, y_center, zoom, C);
    free(map);
    return 0;
}