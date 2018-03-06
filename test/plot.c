#include <stdio.h>

void plot(char* file, double* data, int dim, int size) {
  FILE *gnuplot = popen("/usr/local/bin/gnuplot", "w");
  fprintf(gnuplot, "set terminal png size 800,600\n");
  fprintf(gnuplot, "set output '%s'\n", file);
  fprintf(gnuplot, "set style line 1 lc rgb 'blue' pt 20\n");
  fprintf(gnuplot, "plot '-' w p ls 1\n");
  for (int i = 0; i < size; i++) {
    fprintf(gnuplot, "%e %e\n", data[i * dim], data[i * dim + 1]);
  }
  fprintf(gnuplot, "e\n");
  fflush(gnuplot);
}