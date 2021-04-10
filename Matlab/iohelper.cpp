#include <stdio.h>
#include <math.h>
#include <iostream>

//Defined function round but there is a function round in cmath
//Changed function name to roundd
int roundd(double x){
  return (x - floor(x) > 0.5) ? floor(x)+1 : floor(x);
}


double *normalizeRepresentationIn(double *x, int width, int height){
  double *out = new double[width*height];
  for(int i=0; i<width; i++){
    for(int j=0; j<height; j++){
        out[j*width+i] = x[i*height+j];
    }
  }
  return out;
}

double *normalizeRepresentationOut(double *x, int width, int height){
  double *out = new double[width*height];
  for(int i=0; i<width; i++) {
    for(int j=0; j<height; j++) {
      out[i*height+j] = x[j*width+i];
    }
  }
  return out;
}

double *resize(double *image, int width, int height, double scale){
  int width1 = roundd((double)width*scale);
  int height1 = roundd((double)height*scale);
  double kx = (double)width/(double)width1;
  double ky = (double)height/(double)height1;
  double *resizedImage = new double[width1*height1];

  double alphax, alphay, alpha, alphaSum;
  for(int x=0; x<width1; x++)
    for(int y=0; y<height1; y++){
      resizedImage[y*width1+x]=0;
      alphaSum=0;
      //double alphaSum=0;
      for(int i=floor(x*kx - kx/2); i<ceil(x*kx+kx/2); i++)
        for(int j=floor(y*ky - ky/2); j<ceil(y*ky+ky/2); j++){
          if(i < 0 || i >= width || j < 0 || j >= height)
            continue;
          // double alphax = (i < ceil(x*kx - kx/2)) ? 1 - ((x*kx - kx/2) - floor(x*kx - kx/2)) : (i>floor(x*kx + kx/2)) ? x*kx + kx/2 - floor(x*kx + kx/2) : 1;
          // double alphay = (i < ceil(y*ky - ky/2)) ? 1 - ((y*ky - ky/2) - floor(y*ky - ky/2)) : (i>floor(y*ky + ky/2)) ? y*ky + ky/2 - floor(y*ky + ky/2) : 1;
          // double alpha = alphax*alphay;
          alphax = (i < ceil(x*kx - kx/2)) ? 1 - ((x*kx - kx/2) - floor(x*kx - kx/2)) : (i>floor(x*kx + kx/2)) ? x*kx + kx/2 - floor(x*kx + kx/2) : 1;
          alphay = (i < ceil(y*ky - ky/2)) ? 1 - ((y*ky - ky/2) - floor(y*ky - ky/2)) : (i>floor(y*ky + ky/2)) ? y*ky + ky/2 - floor(y*ky + ky/2) : 1;
          alpha = alphax*alphay;
          resizedImage[y*width1+x] += alpha*image[j*width+i];
          alphaSum += alpha;
        }
      resizedImage[y*width1+x] = (alphaSum==0) ? image[(int)floor(y*ky - ky/2)*width + (int)floor(x*kx - kx/2)] : resizedImage[y*width1+x] / alphaSum;
    }
  return resizedImage;
}
