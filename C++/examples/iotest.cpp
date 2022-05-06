/*  This is some example code using "tifimgio.cpp"
 *
 *  Compile with:
 *  g++ -l tiff -lm -w iotest.cpp -o iotest
 *
 *  Author Jamie (jamie.k.eckstein@gmail.com)
 */

#include <string>
#include <math.h>
#include "../qtiffio.cpp"

int main() {
  char fname_in[] = "../Pics/1hr2741_1.tif";
  char fname_out[] = "../Pics/out.tif";
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
  tifio.write(fname_out, image_scaled, false);

  delete[] tifimg;
  delete[] image_scaled;
  delete[] image_out;
}
