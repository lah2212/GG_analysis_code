/*  This is some example code using "tifimgio.cpp"
 */

#include <string>
#include <math.h>
#include "tifimgio.cpp"

int main() {
  char fname_in[] = "Pics/Pt_94kx_Conical_20um_100_16bit.tif";
  uint32 width, height;

  QTiffIO tifio;
  uint16 *tifimg = tifio.open("Pics/Pt_94kx_Conical_20um_100_16bit.tif");
  tifio.get_dimensions(&width, &height);

  uint16 max = tifimg[0];
  uint16 min = tifimg[0];
  
  for (int i = 0; i < (width * height); i++) {
    if (max < tifimg[i])
      max = tifimg[i];

    if (min > tifimg[i])
      min = tifimg[i];
  }

  double *image_scaled = new double[width * height];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      image_scaled[(width * i) + j] = 65535.0/((double) max - (double) min) *
        (double) (tifimg[(width * i) + j] - min);
    }
  }

  uint16 *image_out = new uint16[width * height];
  for (int i = 0; i < width * height; i++) {
    image_out[i] = (uint16) image_scaled[i];
  }
  
  //QTiffIO tifout;
  tifio.set_dimension(width, height);
  tifio.write("Pics/out.tif", image_scaled, false);

  delete[] tifimg;
  delete[] image_scaled;
  delete[] image_out;
}
