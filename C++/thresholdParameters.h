#include <math.h>


class thresholdParameters{
 private:
  double mean, std;
 public:
  thresholdParameters(double *edgeInfo, int width, int height) {
    mean=0;
    std=0;
    for(int i=0; i<width*height; i++)
      mean += edgeInfo[i];
    mean = mean/(double)(width*height);
    for(int i=0; i<width*height; i++)
      std += pow(edgeInfo[i] - mean, 2);
    std = sqrt(std/(double)(width*height));
  }
  double lowThreshold(){
    return 2.44*std + 0.99*mean - 0.96;
  }
  
  double highThreshold(){
    return 3.6*std + 1.4*(1-mean) - 0.213;
  }
};
