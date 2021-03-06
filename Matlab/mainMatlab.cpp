
//Fixed Error: In function 'int round(double)'
//Also Error in declaration of MexFunction

#include <stdio.h>
#include <iostream>
#include "Smoothing.h"
#include "Morphology.h"
#include "thresholdParameters.h"
#include "interpolation.h"
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
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){ 
  printf("Using %d Images\n", nrhs);
  int width = mxGetN(prhs[0]);
  int height = mxGetM(prhs[0]);  
  int numImages = nrhs;
  vector<double *> imageSet;
//  double scale = 1;
  double scale = 700/sqrt((double)width*height);
  printf("Scale: %f \t Width: %d \t Height: %d\n", scale, int(scale*width), int(scale*height));
  for(int i=0; i<numImages; i++){
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
  delete[] edges;

  Morphology *m = new Morphology();
  double *thresholdedImage = m->doubleThreshold(tmp, tps->lowThreshold(), tps->highThreshold(), width, height);
  delete tps;
  delete[] tmp;

  printf("Denoising Image...\n");
  int firstDenoiseThresh = roundd(width*height/5000);
  printf("Binary Denoising Image...\n");
  double *denoisedImage = m->binaryDenoise(thresholdedImage, width, height, firstDenoiseThresh, 2);
  delete[] thresholdedImage;

  printf("Dilating Image...\n");
  double *dilatedImage = m->dilate(denoisedImage, width, height);
  delete[] denoisedImage;
  printf("Thinning Image...\n");
  double *thinnedImage = m->thin(dilatedImage, width, height);
  delete[] dilatedImage;
//  printf("Detangling Image...\n");
  double *detangledImage = m->detangle(m->detangle(thinnedImage, width, height), width, height);
  delete[] thinnedImage;

  Interpolation *interp = new Interpolation();
  printf("Interpolating Image...\n");
  double *interpolatedImage = interp->interpolate(detangledImage, width, height);
  delete[] detangledImage;
  int pruneThreshold = roundd(sqrt((double)(width*height))/50.0);
  printf("Pruning Image...\n");
  double *prunedImage = m->prune(interpolatedImage, pruneThreshold, width, height);
  delete[] interpolatedImage;
  int finalDenoiseThresh = roundd(sqrt((double)width*height)/50.0);

  double *skeleton = prunedImage;
  double *skeletonDenoised;
  for(int i=0; i<3; i++) {
    printf("Skeletonizing Image...Pass %d of 3\n", i + 1);
    skeletonDenoised = m->binaryDenoise(skeleton, width, height, finalDenoiseThresh, 2);
    delete[] skeleton;
    skeleton = m->prune(m->detangle(interp->interpolate(skeletonDenoised, width, height), width, height), pruneThreshold, width, height);
    delete[] skeletonDenoised;
  }
  printf("Normalizing Image...\n");
  double *outImage = normalizeRepresentationOut(skeleton, width, height);
//  delete[] prunedImage;
  delete[] skeleton;
  delete interp;
  delete m;

  plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
  double *out = mxGetPr(plhs[0]);

  for(int i=0; i<width*height; i++)
    out[i] = outImage[i];

  delete[] outImage;
}
