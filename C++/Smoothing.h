#include <math.h>
#include <stdlib.h>
#include <vector>
#include "qtiffio.cpp"
using namespace std;


class NonlinearIso{
  private:
    double dx(double *image, int x, int y, int width, int height){
      y = (y < 0) ? 0 : ((y > height - 1) ? height - 1 : y);
      x = (x <= 0) ? 1 : ((x >= width - 1) ? width - 2 : x);
      return 0.5 * (image[y * width + x + 1] - image[y * width + x - 1]);
    }
    double dy(double *image, int x, int y, int width, int height){
      x = (x < 0) ? 0 : ((x > width - 1) ? width - 1 : x);
      y = (y <= 0) ? 1 : ((y >= height - 1) ? height - 2 : y);
      return 0.5 * (image[(y + 1) * width + x] - image[(y - 1) * width + x]);
    }
    double dxx(double *image, int x, int y, int width, int height){
      y = (y < 0) ? 0 : ((y > height - 1) ? height - 1 : y);
      x = (x <= 0) ? 1 : ((x >= width - 1) ? width - 2 : x);
      return image[y * width + x + 1] - 2 * image[y * width + x] + image[y * width + x - 1];
    }
    double dyy(double *image, int x, int y, int width, int height){
      x = (x < 0) ? 0 : ((x > width - 1) ? width - 1 : x);
      y = (y <= 0) ? 1 : ((y >= height - 1) ? height - 2 : y);
      return image[(y + 1) * width + x] - 2 * image[y * width + x] + image[(y - 1) * width + x];
    }
    

    double *diffuse(double *image, int width, int height){
      double dt = 0.1;
      double iterations = 30;
      double *diffusedImage = new double[width*height];
      double *tmp = new double[width*height];
      double xx, yy;

      for(int i=0; i<width*height; i++){
        tmp[i] = image[i];
        diffusedImage[i] = image[i];
      }
      
      for(int t=0; t<iterations; t++){
        for(int i=0; i<width; i++)
          for(int j=0; j<height; j++){
            xx=dxx(tmp, i, j, width, height);
            yy=dyy(tmp, i, j, width, height);
            diffusedImage[j*width+i] = tmp[j*width+i] + dt*(xx+yy);
          }
        if(t!=iterations-1)
          for(int i=0; i<width*height; i++)
            tmp[i] = diffusedImage[i];
      }
      delete[] tmp;
      return diffusedImage;
    }

    

    double *transpose(double *x, int width, int height){
      double *out = new double[width*height];
      for(int i=0; i<width; i++)
        for(int j=0; j<height; j++)
            out[i*height+j] = x[j*width+i];
      return out;
    }
    
    
    double *thomas(double *a, double *b, double *c, double *d, int width, int height){
      double *m = new double[width*height];
      double *l = new double[width*height];
      double *y = new double[width*height];

      for(int i=0; i<width*height; i++){
        m[i]=0;
        l[i]=0;
        y[i]=0;
      }

      for(int i=0; i<width; i++) {
        m[i] = a[i];
        y[i] = d[i];
      }

      for(int j=1; j<height; j++) {
        for(int i=0; i<width; i++) {
          l[(j-1)*width+i] = c[(j-1)*width+i]/m[(j-1)*width+i];
          m[j*width+i] = a[j*width+i] - l[(j-1)*width+i]*b[(j-1)*width+i];
          y[j*width+i] = d[j*width+i] - l[(j-1)*width+i]*y[(j-1)*width+i];
        }
      }

      double *x = new double[width*height];
      for(int i=0; i<width; i++)
        x[(height-1)*width+i] = y[(height-1)*width+i]/m[(height-1)*width+i];
      for(int j=height-2; j>=0; j--)
        for(int i=0; i<width; i++)
          x[j*width+i] = (y[j*width+i] - b[j*width+i]*x[(j+1)*width+i])/m[j*width+i];
      
      delete[] m;
      delete[] l;
      delete[] y;

      return x;
    }
    
    
    void *aosiso(double *image, double *tensor, double dt, int width, int height, double *out){
      double *y = new double[width*height];
  //    double *yR = new double[width*height];
  //    double *yC = new double[width*height];
      double *p = new double[width*height];
      double *q = new double[width*height];
      double *a = new double[width*height];
      double *b = new double[width*height];

      for(int i=0; i<width*height; i++) {
        y[i]=0;
  //      yR[i]=0;
  //      yC[i]=0;
        p[i]=0;
        q[i]=0;
        a[i]=0;
        b[i]=0;
      }

      for(int i=0; i<width; i++)
        for(int j=0; j<height-1; j++)
          q[j*width+i] = tensor[j*width+i] + tensor[(j+1)*width+i];

      for(int i=0; i<width; i++) {
        p[i] = q[i];
        p[(height-1)*width+i] = q[(height-2)*width+i];
      }

      for(int i=0; i<width; i++)
        for(int j=1; j<height-1; j++)
          p[j*width+i] = q[(j-1)*width+i] + q[j*width+i];

      for(int i=0; i<width*height; i++){
        a[i] = 1 + dt * p[i];
        b[i] = (-1 * dt) * q[i];
      }
      double *yR = thomas(a, b, b, image, width, height);
      
      for(int i=0; i<width*height; i++){
        p[i]=0;
        q[i]=0;
        a[i]=0;
        b[i]=0;
      }

      for(int j=0; j<height; j++)
        for(int i=0; i<width-1; i++)
          q[j*width+i] = tensor[j*width+i]+tensor[j*width+i+1];

      for(int i=0; i<height; i++){
        p[i*width]=q[i*width];
        p[i*width+width-1] = q[i*width+width-2];
      }

      for(int i=1; i<width-1; i++)
        for(int j=0; j<height; j++)
          p[j*width+i] = q[j*width+i-1]+q[j*width+i];

      double *pT = transpose(p, width, height);
      double *qT = transpose(q, width, height);

      double *xT = transpose(image, width, height);

      for(int i=0; i<width*height; i++){
        a[i] = 1 + dt * pT[i];
        b[i] = (-1*dt) * qT[i];
      }

      double *yT = thomas(a, b, b, xT, height, width);
      double *yC = transpose(yT, height, width);

      for(int i=0; i<width*height; i++)
        y[i]=(yR[i]+yC[i])/2;

      for (int i = 0; i < width * height; i++) {
        out[i] = y[i];
      }

      delete[] a;
      delete[] b;
      delete[] pT;
      delete[] qT;
      delete[] p;
      delete[] q;
      delete[] yT;
      delete[] yC;
      delete[] yR;
      delete[] xT;
      
      delete[] y;
  //    return y;
    }

    /*
    double LoGFilter[] =  // sigma = 1.4
    {0.0002, 0.0008, 0.0021, 0.0036, 0.0043, 0.0036, 0.0021, 0.0008, 0.0002, 
     0.0008, 0.0030, 0.0070, 0.0100, 0.0108, 0.0100, 0.0070, 0.0030, 0.0008, 
     0.0021, 0.0070, 0.0112, 0.0064, 0.0006, 0.0064, 0.0112, 0.0070, 0.0021, 
     0.0036, 0.0100, 0.0064, -0.0244, -0.0478, -0.0244, 0.0064, 0.0100, 0.0036, 
     0.0043, 0.0108, 0.0006, -0.0478, -0.0829, -0.0478, 0.0006, 0.0108, 0.0043, 
     0.0036, 0.0100, 0.0064, -0.0244, -0.0478, -0.0244, 0.0064, 0.0100, 0.0036, 
     0.0021, 0.0070, 0.0112, 0.0064, 0.0006, 0.0064, 0.0112, 0.0070, 0.0021, 
     0.0008, 0.0030, 0.0070, 0.0100, 0.0108, 0.0100, 0.0070, 0.0030, 0.0008, 
     0.0002, 0.0008, 0.0021, 0.0036, 0.0043, 0.0036, 0.0021, 0.0008, 0.0002}
    */

    double *lapFilter(double *image, int width, int height, bool freeimage){
      double laplacianFilter[] = {0.16667, 0.66667, 0.16667,
                                  0.66667, -3.3333, 0.66667,
                                  0.16667, 0.66667, 0.16667};

      double *output = filter2(image, width, height, laplacianFilter, 3);
      if (freeimage) {
          delete[] image;
      }
      return output;
    }
    
    double *medFilter(double *image, int width, int height, bool freeimage){
      double medianFilter[] = {0.11111, 0.11111, 0.11111, 
                               0.11111, 0.11111, 0.11111, 
                               0.11111, 0.11111, 0.11111}; 

      double *output = filter2(image, width, height, medianFilter, 3);
      if (freeimage) {
          delete[] image;
      }
      return output;
    }
    
    
    // So I'm pretty sure this is just convolving
    double *filter2(double *image, int width, int height, double* filter, int fWidth){
      double *filteredImage = new double[width*height];
      int q = floor((double)(fWidth)/2.0);

      for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
          filteredImage[j*width+i] = 0;
          for(int m=i-q; m<=i+q; m++) {
            for(int n=j-q; n<=j+q; n++) {
              int x = (m<0) ? -1*m : ((m>=width) ? width - (1 + m - width) : m);
              int y = (n<0) ? -1*n : ((n>=height) ? height - (1 + n - height) : n);
              int xx = m - i + q;
              int yy = n - j + q;
              filteredImage[j*width+i] += image[y*width+x] * filter[yy*fWidth+xx];
            }
          }
        }
      }

      return filteredImage;
    }
   
   
    
    double *copy(double *tmp, int width, int height){
      double *out = new double[width*height];
      for(int i=0; i<width*height; i++)
        out[i] = tmp[i];
      return out;
    }
   
    double *diffusionTensor(vector<double *> images, int width, int height, double k, bool filterFirst, bool verbose = false) {
      double *tensor = new double[width*height];
      for(int i = 0; i < width * height; i++)
        tensor[i] = 0;

      double *medFilteredImage;
      for(int imageNum=0; imageNum<images.size(); imageNum++) {
        if(filterFirst)
          medFilteredImage = medFilter(images[imageNum], width, height, false);
        else {
          medFilteredImage = new double[width*height];
          for(int i=0; i<width*height; i++)
            medFilteredImage[i] = images[imageNum][i];
  //         delete[] medFilteredImage;
  //         double *medFilteredImage = copy(images[imageNum], width, height);
        }

        for(int i=0; i<width; i++)
          for(int j=0; j<height; j++)
            tensor[j*width+i] += pow(sqrt(pow(dx(medFilteredImage, i, j, width, height), 2) + 
                              pow(dy(medFilteredImage, i, j, width, height), 2)), 2);
      }
      delete[] medFilteredImage;

      if (!verbose) {
        for(int i=0; i<width*height; i++)
          tensor[i] = exp(-k*sqrt(tensor[i]));

        tensor = medFilter(tensor, width, height, true);

        tensor = lapFilter(tensor, width, height, true);

        double *tmp = new double[width*height];
        for(int i=0; i<width*height; i++)
          tmp[i] = (tensor[i]>0) ? tensor[i] : 0;
        double *tmpMed = medFilter(tmp, width, height, false);
        for(int i=0; i<width*height; i++)
          tensor[i] = exp(-10*tmpMed[i]);

        delete[] tmp;
        delete[] tmpMed;
      }
      else {
        std::printf("Verbose\n");
        QTiffIO tensor_tif;
        tensor_tif.set_dimension(width, height);

        tensor_tif.write("Pics/tensor0.tif", tensor, true);
        
        for(int i=0; i<width*height; i++)
          tensor[i] = exp(-k*sqrt(tensor[i]));

        tensor_tif.write("Pics/tensor2.tif", tensor, true);
        
        tensor = medFilter(tensor, width, height, true);

        tensor_tif.write("Pics/tensor3.tif", tensor, true);
        
        tensor = lapFilter(tensor, width, height, true);

        tensor_tif.write("Pics/tensor4.tif", tensor, true);

        double *tmp = new double[width*height];
        for(int i=0; i<width*height; i++)
          tmp[i] = (tensor[i]>0) ? tensor[i] : 0;
        double *tmpMed = medFilter(tmp, width, height, false);
        for(int i=0; i<width*height; i++)
          tensor[i] = exp(-10*tmpMed[i]);

        tensor_tif.write("Pics/tensor5.tif", tensor, true);

        delete[] tmp;
        delete[] tmpMed;
      }

      return tensor;
    }



   public:
     
    double *nonlinearIso(vector<double *> images, double dt, int iterations, int width, int height) {
      vector<double *> diffImages;
      int numImages = images.size();

      for(int i=0; i<numImages; i++)
        diffImages.push_back(medFilter(images[i], width, height, false));
        
      char buffer[80];
      for(int t = 0; t < iterations; t++){
        double *diffTens = diffusionTensor(diffImages, width, height, 0.08+t/(12*iterations), (t%2==0)?true:false);
        for(int k = 0; k < numImages; k++) {
          //diffImages[k] = aosiso(images[k], diffTens, dt, width, height, diffImages[k]);
          aosiso(images[k], diffTens, dt, width, height, diffImages[k]);
        }
        sprintf(buffer, "Pics/difftens_iter_%d.tif", t);
        delete[] diffTens;
      }

    double *diffTens = diffusionTensor(diffImages, width, height, 0.1633, false, true);

    /*
    diffTens = medFilter(diffTens, width, height, true);
    diffTens = lapFilter(diffTens, width, height, true);
    double *tmp = new double[width*height];
    for(int i=0; i<width*height; i++)
        tmp[i] = (diffTens[i]>0) ? diffTens[i] : 0;
    double *tmpMed = medFilter(tmp, width, height, false);
    for(int i=0; i<width*height; i++)
        diffTens[i] = exp(-10*tmpMed[i]);
    
    delete[] tmp;
    delete[] tmpMed;
      */

    for(int i = 0; i < numImages; i++)
      delete[] diffImages[i];

    return diffTens;
  }
  
};
