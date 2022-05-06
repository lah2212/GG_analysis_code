/*  Feeds n-channel image stack into mainMatlab.cpp algorithm
 *
 *  Author:  Jamie Eckstein (), jamie.k.eckstein@gmail.com
 */
#include <stdlib.h>
#include <iostream>
#include <string>
//#include "qtiffio.cpp"
#include "mainMatlab.cpp"

int main(int argc, char **argv) {
  const int NUM_IMAGES = 3;
//  const double SCALE = 700.0/2048.0;
//  const double SCALE = 1;
  const double SCALE = 700/sqrt(4012.0 * 4012.0);
  
  double *image_stack[NUM_IMAGES];
  QTiffIO tifio;
  uint32 width, height;

//  char* fnames[1] = { "Pics/STEM/Pt170_STEM_225kX_C2(100)_CL205mm_06.tif" };

  /*
  char* fnames[4] = {
        "Pics/1hr2741_1.tif",
        "Pics/1hr2741_2.tif",
        "Pics/1hr2741_3.tif",
        "Pics/1hr2741_4.tif"};
  char *fnames[4] = {
        "Pics/1hr2789_1.tif",
        "Pics/1hr2789_2.tif",
        "Pics/1hr2789_3.tif",
        "Pics/1hr2789_4.tif"};
  */
  char *fnames[3] = {
        "Pics/2hr2347_1.tif",
        "Pics/2hr2348_2.tif",
        "Pics/2hr2349_3.tif"};

  for(int i = 0; i < NUM_IMAGES; i++) {
    std::printf("Opening %s\n", fnames[i]);
    uint16 *image = tifio.open(fnames[i]);
    tifio.get_dimensions(&width, &height);

    uint16 max = image[0];
    uint16 min = image[0];
    
    for (int i = 0; i < (width * height); i++) {
      if (max < image[i])
        max = image[i];

      if (min > image[i])
        min = image[i];
    }

//    std::printf("Max: %d, Min: %d\n", max, min);

    double *image_scaled = new double[width * height];
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        image_scaled[(width * i) + j] = 255.0/((double) max - (double) min) *
          (double) (image[(width * i) + j] - min);
      }
    }
    delete[] image;

    image_stack[i] = image_scaled;
  }

  double *image_processed = 
               algorithm(NUM_IMAGES, image_stack, width, height, SCALE);

  int width_out = SCALE * width;
  int height_out = SCALE * height;

  tifio.set_dimension(width_out, height_out);
  tifio.write("Pics/finalskel.tif", image_processed, true);

  return 0;
}
