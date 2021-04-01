/*
 *       Filename:  main.cpp
 *    Description:  Just the file that handles all the image I/O
 *        Version:  1.0
 *         Author:  Jamie Eckstein (), jamie.k.eckstein@gmail.com
 */
#include <stdlib.h>
#include <iostream>
#include <string>
//#include "qtiffio.cpp"
#include "mainMatlab.cpp"

int main(int argc, char **argv) {
  const int NUM_IMAGES = 5;
//  const double SCALE = 700.0/2048.0;
  const double SCALE = 1;
//  const double SCALE = 700/sqrt(4096.0 * 4096.0);
  
  double *image_stack[NUM_IMAGES];
  QTiffIO tifio;
  uint32 width, height;

//  char* fnames[1] = { "Pics/STEM/Pt170_STEM_225kX_C2(100)_CL205mm_06.tif" };

  char* fnames[5] = {
        "Pics/TEM/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_96_16bit.tif",
        "Pics/TEM/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97_16bit.tif",
        "Pics/TEM/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_98_16bit.tif",
        "Pics/TEM/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_99_16bit.tif",
        "Pics/TEM/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_100_16bit.tif"};

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
  tifio.write("finalskel.tif", image_processed, true);

  return 0;
}
