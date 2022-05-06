/*  This just handles the image io into mainMatlab.cpp using qtiffio
 *
 *  Author:  Jamie Eckstein, (jamie.k.eckstein@gmail.com)
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include "mainMatlab.cpp"

int main(int argc, char **argv) {
  char *fname_in = "Pics/10hr2144_1.tif";
  if (argc <= 1) {
    std::printf("%s\n", fname_in);
  } 
  else {
    fname_in = argv[1];
    std::printf("%s\n", fname_in);
  }

  char fname_out[] = "Pics/out.tif";
  uint32 width, height;

  QTiffIO tifio;
  uint16 *tifimg = tifio.open(fname_in, true);
  tifio.get_dimensions(&width, &height);

  uint16 max = tifimg[0];
  uint16 min = tifimg[0];
  
  for (int i = 0; i < (width * height); i++) {
    if (max < tifimg[i])
      max = tifimg[i];

    if (min > tifimg[i])
      min = tifimg[i];
  }

  std::printf("Max: %d, Min: %d\n", max, min);

  double *image_scaled = new double[width * height];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      image_scaled[(width * i) + j] = 65535.0/((double) max - (double) min) *
        (double) (tifimg[(width * i) + j] - min);
    }
  }

  double *image_stack[1];
  image_stack[0] = image_scaled;
  double scale = 700/sqrt((double)width*height);
  double *image_processed = algorithm(1, image_stack, width, height, scale);

  int width_out = scale * width;
  int height_out = scale * height;

  double max_s = image_processed[0];
  double min_s = image_processed[0];
  
  for (int i = 0; i < (width_out * height_out); i++) {
    if (max_s < image_processed[i])
      max_s = image_processed[i];

    if (min_s > image_processed[i])
      min_s = image_processed[i];
  }
  std::printf("Max: %f, Min: %f\n", max_s, min_s);

  uint16 *image_out = new uint16[width_out * height_out];
  for (int i = 0; i < width_out * height_out; i++) {
    image_out[i] = (uint16) 65535 * image_processed[i];
  }

  tifio.set_dimension(width_out, height_out);
  tifio.write(fname_out, image_out);

  delete[] tifimg;
  delete[] image_scaled;
  delete[] image_processed;
  delete[] image_out;
  return 0;
}
