#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <sys/mman.h>

unsigned char *mandelbrot_mapping_mpi(size_t N, int portion, int begin, int end, size_t C, double inc, double x_min, double y_max)
{
    double limit = 2.0;
    unsigned char *g_map = aligned_alloc(64, portion * N* sizeof(unsigned char));
    unsigned char iteration = (unsigned char)C;
    unsigned char cutoff = 0;
    memset(g_map, 0, portion * N * sizeof(unsigned char));

    for (int y = begin; y < end; y++)
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
            g_map[(y-begin) * N + x] = iteration -cutoff;
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
  // First initialize MPI
  MPI_Init(&argc , &argv);

  // Initialize variables for communicators
  int comm_size, comm_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  

  size_t N, C; // order, cutoff
  double x_center, y_center, zoom;
  if (argc < 6)
  {
    fprintf(stderr,
            "Usage: ./mandelbrot_serial order xcenter ycenter zoom cutoff\n");
    exit(1);
  }
  sscanf(argv[1], " %zu ", &N);
  sscanf(argv[2], " %lf ", &x_center);
  sscanf(argv[3], " %lf ", &y_center);
  sscanf(argv[4], " %lf ", &zoom);
  sscanf(argv[5], " %zu ", &C);
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

  double inc = pow(2, (-1) * zoom);
  double width = inc * N;
  double x_min = x_center - width / 2;
  double y_max = y_center + width / 2;

  int rec_array[2];
  int portion = N / comm_size;
  unsigned char *final_map = NULL;
  int * task_array = NULL;

  if(!comm_rank) { // operations only for p0
      task_array = malloc(comm_size * 2 * sizeof(int));
      for(int i = 0; i < comm_size*2; i++) {
          task_array[i] = (i/2) * portion;
          task_array[i+1] = ((i/2)+1) * portion;
          i++;
      }
      final_map = aligned_alloc(64, N * N * sizeof(unsigned char));
  }
  
  // Every process receives an array of size 2, that contains the beginning and ending index of their task
  MPI_Scatter(task_array , 2 , MPI_UNSIGNED , rec_array , 2 , MPI_UNSIGNED , 0 , MPI_COMM_WORLD);

  int begin = rec_array[0];
  int end = rec_array[1];
  unsigned char *local_map = mandelbrot_mapping_mpi(N, portion, begin, end, C, inc, x_min, y_max);

  MPI_Gather(local_map , portion * N , MPI_UNSIGNED_CHAR , final_map , portion * N , MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if(!comm_rank) {
    write_image(final_map, N, x_center, y_center, zoom, C);
  }

  MPI_Finalize();
  return 0;
}