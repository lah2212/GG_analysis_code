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

int main(int argc, char **argv) {
    char *pic = new char[100];
    if(argc <= 1) {
      pic = "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif";
    } else {
      pic = argv[1];
    }

    std::printf("%s\n", pic);

    TIFF *tif = TIFFOpen(pic, "r");
    if(tif) {
      std::printf("Toite\n");

      TIFF *out = TIFFOpen("Pics/skel.tif", "w");
//        out = tif;
      TIFFClose(out);
    }
    TIFFClose(tif);

    return 0;
}
