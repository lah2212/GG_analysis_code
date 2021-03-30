/*  This is just a little program for importing and exporting tif files.
 *  Designed to only work with the PT_Conical images downsized to 16 bit.
 *
 *  Author:  Jamie Eckstein (), jamie.k.eckstein@gmail.com
 */

#include <iostream>
#include <string>
#include "tiffio.h"

class QTiffIO {
    uint32 width, height;
    uint32 width_out, height_out;
    bool fetched;

  public:
    QTiffIO() {
      fetched = false;
    }

    void get_dimensions(uint32 *set_width, uint32 *set_height) {
      if (fetched) {
        *set_width = width;
        *set_height = height;
      }
      else {
        std::printf("Image needs to be opened before dimensions can be gathered\n");
        exit(1);
      }
    }

    uint16* open(const char *path, bool verbose = false) {
      uint16 *image;
      TIFF *tif;
      tif = TIFFOpen(path, "r");

      if (tif) {
        uint16 nsamples;
        uint32 config;
        fetched = true;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
        if (verbose) {
          std::printf("Opening Tiff File...\n");
          std::printf("Samples: %d\n", nsamples);
          std::printf("Image (width, height): (%d, %d)\n", width, height);
        }

        image = new uint16[width * height];
        uint16 *buf = (uint16 *) _TIFFmalloc(TIFFScanlineSize(tif));
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
        return image;
      } 
      else {
        std::printf("Could not open tiff image");
        exit(1);
      }
    }

    void set_dimension (uint32 width, uint32 height) {
      width_out = width;
      height_out = height;
    }

    void write(const char *path, double *image, bool normalize) {
      uint16 *image_out = new uint16[width_out * height_out];

      if (normalize) {
        double max = image[0];
        double min = image[0];
        
        for (int i = 0; i < (width_out * height_out); i++) {
          if (max < image[i])
            max = image[i];

          if (min > image[i])
            min = image[i];
        }

        std::printf("Max: %f, Min: %f\n", max, min);
        double *image_scaled = new double[width_out * height_out];
        for (int i = 0; i < height_out; i++) {
          for (int j = 0; j < width_out; j++) {
            image_scaled[(width_out * i) + j] = 65535.0/((double) max - (double) min) *
              (double) (image[(width_out * i) + j] - min);
          }
        }

        for (int i = 0; i < width_out * height_out; i++) {
          image_out[i] = (uint16) image_scaled[i];
        }

        delete[] image_scaled;
      }
      else {
        for (int i = 0; i < width_out * height_out; i++) {
          image_out[i] = (uint16) image[i];
        }
      }

      write(path, image_out);
      delete[] image_out;
    }

    void write(const char *path, uint16 *image_out) {
      std::printf("Writing image...\n");
      TIFF *out = TIFFOpen(path, "w");
      if (out) {
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width_out);  // set the width of the image
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, height_out);    // set the height of the image
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 16);    // set the size of the channels
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

        TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, height));
        TIFFSetupStrips(out); 

        tsize_t linebytes = width_out * 2;

        uint16 *buf;
        if (TIFFScanlineSize(out) == linebytes) {
          buf = (uint16 *) _TIFFmalloc(linebytes);
        }
        else {
          buf = (uint16 *) _TIFFmalloc(TIFFScanlineSize(out));
        }
        
        //Now writing image to the file one strip at a time
        for (uint32 row = 0; row < height_out; row++) {
  //        memcpy(buf, &image_out[(height_out - row - 1) * linebytes], linebytes);    // check the index here, and figure out why not using h*linebytes
          memcpy(buf, &image_out[(row) * width_out], linebytes);    // check the index here, and figure out why not using h*linebytes
          if (TIFFWriteScanline(out, buf, row, 0) < 0)
            break;
        }

        std::printf("Image Closed...\n");
        (void) TIFFClose(out);
        _TIFFfree(buf);
      }
      else {
        std::printf("Failed to open tif %s for exporting\n", path);
        exit(1);
      }
    }
};
