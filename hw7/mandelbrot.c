//Branden Foley, bpfoley, modified mandlebrot.c, parallelize(?) this process with openMP.


// Turn in stuff
//   1. 4 cores (lab machine)
//   2. times




/*
  This program is an adaptation of the Mandelbrot program
  from the Programming Rosetta Stone, see
  http://rosettacode.org/wiki/Mandelbrot_set

  Compile the program with:

  gcc -o mandelbrot -O4 mandelbrot.c

  Usage:

  ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>

  Example:

  ./mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm

  The interior of Mandelbrot set is black, the levels are gray.
  If you have very many levels, the picture is likely going to be quite
  dark. You can postprocess it to fix the palette. For instance,
  with ImageMagick you can do (assuming the picture was saved to pic.ppm):

  convert -normalize pic.ppm pic.png

  The resulting pic.png is still gray, but the levels will be nicer. You
  can also add colors, for instance:

  convert -negate -normalize -fill blue -tint 100 pic.ppm pic.png

  See http://www.imagemagick.org/Usage/color_mods/ for what ImageMagick
  can do. It can do a lot.
*/

//TODO --- store all the stuff to write to the file, and then only write it after all of the parellizing takes place.


// This code consists of one major loop nest that can be parallelized, but it has a problem. Every iteration is writing into the file, and it must be written in order. You will need to save all the entries somehow and then write them all out after the computation. I strongly suggest you get this working before you parallelize the code.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

int main(int argc, char* argv[])
{
  /* Parse the command line arguments. */
  if (argc != 8) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* The window in the plane. */
  const double xmin = atof(argv[1]);
  const double xmax = atof(argv[2]);
  const double ymin = atof(argv[3]);
  const double ymax = atof(argv[4]);

  /* Maximum number of iterations, at most 65535. */
  const uint16_t maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  const int xres = atoi(argv[6]);
  const int yres = (xres*(ymax-ymin))/(xmax-xmin);

  /* The output file name */
  const char* filename = argv[7];

  /* Open the file and write the header. */
  FILE * fp = fopen(filename,"wb");
  char *comment="# Mandelbrot set";/* comment should start with # */

  /*write ASCII header to the file*/
  fprintf(fp,
          "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
          xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));


//Start stuff here store the  pixels in a big array, then write them all in a proper order.
  typedef unsigned char Pixel[6];
  Pixel *pixels = malloc(sizeof(Pixel) *xres*yres); //this will maybe work?? was different on slack.

  /* Precompute pixel width and height. */
  double dx=(xmax-xmin)/xres;
  double dy=(ymax-ymin)/yres;

  clock_t begin = clock();
  #pragma omp parallel shared(pixels)
  {
  double x, y; /* Coordinates of the current point in the complex plane. */
  double u, v; /* Coordinates of the iterated point. */
  int i,j; /* Pixel counters */
  int k; /* Iteration counter */

  #pragma omp for schedule(dynamic)
  for (j = 0; j < yres; j++) {
    y = ymax - j * dy;
    for(i = 0; i < xres; i++) {
      double u = 0.0;
      double v= 0.0;
      double u2 = u * u;
      double v2 = v*v;
      x = xmin + i * dx;
      /* iterate the point */
      for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
            v = 2 * u * v + y;
            u = u2 - v2 + x;
            u2 = u * u;
            v2 = v * v;
      };
      /* compute  pixel color and write it to file */
      if (k >= maxiter) {
        /* interior */ //NOTE OLD CODE
        // const unsigned char black[] = {0, 0, 0, 0, 0, 0};
        // fwrite (black, 6, 1, fp);

        pixels[j*xres + i][0] = 0;
        pixels[j*xres + i][1] = 0;
        pixels[j*xres + i][2] = 0;
        pixels[j*xres + i][3] = 0;
        pixels[j*xres + i][4] = 0;
        pixels[j*xres + i][5] = 0;
      }
      else {
        /* exterior */ //NOTE old code
        // unsigned char color[6];
        // color[0] = k >> 8;
        // color[1] = k & 255;
        // color[2] = k >> 8;
        // color[3] = k & 255;
        // color[4] = k >> 8;
        // color[5] = k & 255;
        // fwrite(color, 6, 1, fp);

        pixels[j*xres + i][0] = k >>8;
        pixels[j*xres + i][1] = k & 255;
        pixels[j*xres + i][2] = k >>8;
        pixels[j*xres + i][3] = k & 255;
        pixels[j*xres + i][4] = k >>8;
        pixels[j*xres + i][5] = k & 255;
      };
    }
  }
}
  clock_t end = clock();
  double time_spent = (double)(end-begin) /CLOCKS_PER_SEC;
  printf("time: %f\n", time_spent);
  //write to file.
  fwrite(pixels, sizeof(Pixel), xres * yres, fp);
  fclose(fp);
  return 0;
}
