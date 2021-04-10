// This function is just mainMatlab.cpp without the edge detection or smoothing
// This is so smoothing algorithm written in matlab can be tested

#include <stdio.h>
#include <iostream>
#include "qtiffio.cpp"
#include "iohelper.cpp"
#include "Morphology.h"
#include "thresholdParameters.h"
#include "interpolation.h"
#include "mex.h"
#include <vector>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) { 
  if (nrhs > 1) {
    printf("Error...Only one edge image can be processed at a time,\n \
      all n-channel processing must take place earlier.");
    exit(1);
  }

  int width = mxGetN(prhs[0]);
  int height = mxGetM(prhs[0]);  

  QTiffIO tifio;
  tifio.set_dimension(width, height);

  printf("Width: %d \n", width);
  printf("Height: %d \n", height);

  double *edges = normalizeRepresentationIn(mxGetPr(prhs[0]), width, height);

  printf("Threshold Parameters...\n");
  thresholdParameters *tps = new thresholdParameters(edges, width, height);
  double *tmp = new double[width*height];
  for(int i = 0; i < width * height; i++)
    tmp[i] = 1-edges[i];
  delete[] edges;

  printf("Low Threshold: %f\n", tps->lowThreshold());
  printf("High Threshold: %f\n", tps->highThreshold());
  Morphology *m = new Morphology();
  double *thresholdedImage = m->doubleThreshold(tmp, tps->lowThreshold(), tps->highThreshold(), width, height);
  delete tps;
  delete[] tmp;

  tifio.write("Pics/1_threshold.tif", thresholdedImage, true);

  printf("Denoising Image...\n");
  int firstDenoiseThresh = roundd(width*height/5000);
  printf("Binary Denoising Image...\n");
  double *denoisedImage = m->binaryDenoise(thresholdedImage, width, height, firstDenoiseThresh, 2);
  delete[] thresholdedImage;

  tifio.write("Pics/2_denoise.tif", denoisedImage, true);

  printf("Dilating Image...\n");
  double *dilatedImage = m->dilate(denoisedImage, width, height);
  delete[] denoisedImage;

  tifio.write("Pics/3_dilated.tif", dilatedImage, true);

  printf("Thinning Image...\n");
  double *thinnedImage = m->thin(dilatedImage, width, height);
  delete[] dilatedImage;
//  printf("Detangling Image...\n");
  double *detangledImage = m->detangle(m->detangle(thinnedImage, width, height), width, height);
  delete[] thinnedImage;

  /* temp */
  double *outImage = normalizeRepresentationOut(detangledImage, width, height);
  plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
  double *out = mxGetPr(plhs[0]);

  for(int i=0; i<width*height; i++)
    out[i] = outImage[i];

  delete[] outImage;

  /*
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
  */
}
