//Fixed Error: In function 'int round(double)'
//Also Error in declaration of MexFunction

#include <stdio.h>
#include <iostream>
#include "../Smoothing.h"
#include "../Morphology.h"
#include "../thresholdParameters.h"
#include "../interpolation.h"
#include "mex.h"
#include <vector>
using namespace std;

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

//Change to MexFunction
//void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[])
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) { 
  int width = mxGetN(prhs[0]);
  int height = mxGetM(prhs[0]);  
  int numImages = nrhs - 1;
  vector<double *> imageSet;

  int nargs = mxGetN(prhs[nrhs - 1]);
  double *args = mxGetPr(prhs[nrhs - 1]);
  double scale;
  if (args[0] != 0) {
    scale = args[0];
  } else if (nargs == 3) {
    // Add threshold stuff here
  } else {
    scale = 1;
//  scale = 700/sqrt((double)width*height);
  }

  printf("Using %d Images\n", numImages);
  printf("Scale: %f \t Width: %d \t Height: %d\n", scale, int(scale*width), int(scale*height));
  for(int i=0; i < numImages; i++){
    double *in = normalizeRepresentationIn(mxGetPr(prhs[i]), width, height);
    imageSet.push_back(resize(in, width, height, scale));
    delete[] in;
  }

  int iterations = 5;
  width = roundd(width*scale);
  height = roundd(height*scale);
  printf("Beginning to smooth...\n");
  NonlinearIso *ni = new NonlinearIso();
  double *edges = ni->nonlinearIso(imageSet, 10, iterations, width, height);

  for(int i = 0; i < numImages; i++) {
    delete[] imageSet.back();
    imageSet.pop_back();
  }

  delete ni;
  printf("Threshold Parameters...\n");
  thresholdParameters *tps = new thresholdParameters(edges, width, height);
  double *tmp = new double[width*height];
  for(int i = 0; i < width * height; i++)
    tmp[i] = 1-edges[i];
//  delete[] edges;

  Morphology *m = new Morphology();
  double *thresholdedImage = m->doubleThreshold(tmp, tps->lowThreshold(), tps->highThreshold(), width, height);
  delete tps;
  delete[] tmp;

  int firstDenoiseThresh = roundd(width*height/5000);
  printf("Binary Denoising Image...\n");
  double *denoisedImage = m->binaryDenoise(thresholdedImage, width, height, firstDenoiseThresh, 2);
//  delete[] thresholdedImage;

  double *edgeoutImage = normalizeRepresentationOut(edges, width, height);
  double *threshoutImage = normalizeRepresentationOut(thresholdedImage, width, height);
  double *denoiseoutImage = normalizeRepresentationOut(denoisedImage, width, height);

  plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(height, width, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(height, width, mxREAL);
  double *edgeout = mxGetPr(plhs[0]);
  double *threshout = mxGetPr(plhs[1]);
  double *denoiseout = mxGetPr(plhs[2]);

  for(int i=0; i<width*height; i++) {
    edgeout[i] = edgeoutImage[i];
    threshout[i] = threshoutImage[i];
    denoiseout[i] = denoiseoutImage[i];
  }

  delete[] edges;
  delete[] thresholdedImage;
  delete[] denoisedImage;
  delete[] edgeoutImage;
  delete[] threshoutImage;
  delete[] denoiseoutImage;
}
