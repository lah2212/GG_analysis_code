/*
 *       Filename:  main.cpp
 *    Description:  Just the file that handles all the image I/O
 *        Version:  1.0
 *         Author:  Jamie Eckstein (), jamie.k.eckstein@gmail.com
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include "tiffio.h"
#include "mainMatlab.cpp"

int main(int argc, char **argv) {
  //char *pic = new char[100];
  //char pic[100];
  
  uint16 *image;
  uint32 width, height;
  TIFF *tif;
  
  if (argc <= 1) {
    std::printf("%s\n", "Pics/20um_97_16bit.tif");
    tif = TIFFOpen("Pics/Pt_94kx_Conical_20um_100_16bit.tif", "r");
  } 
  else {
    std::printf("%s\n", argv[1]);
    tif = TIFFOpen(argv[1], "r");
  }

  if (tif) {
    uint16 nsamples;
    uint32 config;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
    std::printf("%d\n", nsamples);

    image = new uint16[width * height];
    uint32 *buf = (uint32 *) _TIFFmalloc(TIFFScanlineSize(tif));
    uint16 *data;

    for (uint32 row = 0; row < height; row++) {
      if (TIFFReadScanline(tif, buf, row)) {
        data = (uint16 *) buf;
        for (int i = 0; i < width; i++) {
          image[(width * row) + i] = data[i];
        }
      }
      else {
        std::printf("Couldn't read scanlines\n");
        exit(1);
      }
    }

    _TIFFfree(buf);
    TIFFClose(tif);
  } 
  else {
    std::printf("Could not open tiff image");
    exit(1);
  }

  uint16 max = image[0];
  uint16 min = image[0];
  
  for (int i = 0; i < (width * height); i++) {
    if (max < image[i])
      max = image[i];

    if (min > image[i])
      min = image[i];
  }

  std::printf("Max: %d, Min: %d\n", max, min);

  double *image_scaled = new double[width * height];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      image_scaled[(width * i) + j] = 255.0/((double) max - (double) min) *
        (double) (image[(width * i) + j] - min);
    }
  }
  delete[] image;


//    double *image_processed = image_scaled;
  double *image_stack[1];
  image_stack[0] = image_scaled;
  double scale = 700/sqrt((double)width*height);
//  double scale = 1;
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

  uint8 *image_out = new uint8[width_out * height_out];
  for (int i = 0; i < width_out * height_out; i++) {
    image_out[i] = (uint8) 255 * image_processed[i];
  }
  delete[] image_scaled;
  delete[] image_processed;

  TIFF *out = TIFFOpen("new.tif", "w");
  TIFFSetField (out, TIFFTAG_IMAGEWIDTH, width_out);  // set the width of the image
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, height_out);    // set the height of the image
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, height));
  TIFFSetupStrips(out); 

  tsize_t linebytes = width_out;

  unsigned char *buf;
//    uint16 *buf = (uint16 *) _TIFFmalloc(TIFFScanlineSize(out));
  if (TIFFScanlineSize(out) == linebytes)
    buf = (unsigned char *) _TIFFmalloc(linebytes);
  else
    buf = (unsigned char *) _TIFFmalloc(TIFFScanlineSize(out));
  
  //Now writing image to the file one strip at a time
  for (uint32 row = 0; row < height_out; row++) {
    memcpy(buf, &image_out[(height_out - row - 1) * linebytes], linebytes);    // check the index here, and figure out why not using h*linebytes
    if (TIFFWriteScanline(out, buf, row, 0) < 0)
      break;
  }

  (void) TIFFClose(out);
  _TIFFfree(buf);

  delete[] image_out;
  return 0;
}
