#include <stdlib.h>
#ifndef POINT
#define POINT
#include "Point.h"
#endif

#ifndef CORNERDETECTOR
#define CORNERDETECTOR
#include "cornerDetector.h"
#endif

#ifndef ENDDETECTOR
#define ENDDETECTOR
#include "endDetector.h"
#endif

#ifndef TRIPLEDETECTOR
#define TRIPLEDETECTOR
#include "tripleDetector.h"
#endif

#define ptsForSlope 6

#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <queue>
using namespace std;



class LineFunction{
 public:
  double m, b;
  bool bad;
  LineFunction(double m1, double b1){
    m=m1;
    b=b1;
    bad=false;
  }
  LineFunction(bool b){
    bad=b;
  }
  double yGivenX(int x){
    return m*x+b;
  }
  double xGivenY(int y){
    return (y-b)/m;
  }
};

class Morphology{
 private:
  double pi;
  bool elementOfTrip(vector<Point *> v, Point *p, int windowWidth){
    for(int i=0; i<v.size(); i++){
      int x=p->x, y=p->y;
      for(int m=x-windowWidth; m<=x+windowWidth; m++)
	for(int n=y-windowWidth; n<=y+windowWidth; n++){
	  Point *tmp = (Point *)v[i];
	  if(m==tmp->x && n==tmp->y)
	    return true;
	}
    }
    return false;
  }
  
  double myAtanDiff(double x1, double x2){
    double tmp1 = atan(x1);
    double tmp2 = atan(x2);
    return min((double) abs(tmp1-tmp2), 2*pi-abs(tmp1-tmp2));
  }
  
  
  LineFunction *functionOfEndPoint(double *skel, Point *p, int width, int height, vector<Point *> trips){
    vector<Point *> pointList;
    Point *tripPoint=NULL;
    Point *point = new Point(p);
    pointList.push_back(new Point(point));
    while(pointList.size()<ptsForSlope && point!=NULL) {
      int x=point->x;
      int y=point->y;
      point=NULL;
      bool breakLoops=false;
      for(int m=x-1; m<=x+1; m++) {
	for(int n=y-1; n<=y+1; n++) {
	  if((m==x && n==y) || m<0 || m>=width || n<0 || n>=height)
	    continue;
	  if(skel[n*width+m]==0)
	    continue;
	  if(elementOf(pointList, new Point(m, n)))
	    continue;
	  if(elementOfTrip(trips, new Point(m, n), 1)){
	    tripPoint = new Point(m, n);
	    breakLoops=true;
	    break;
	  }
	  point = new Point(m, n);
	  pointList.push_back(new Point(point));
	  breakLoops=true;
	  break;
	}
	if(breakLoops)
	  break;
      }
    }
    int xSum=0, ySum=0, x2=0, xy=0;
    int n=pointList.size();
    int negX=0, posX=0;
    int negY=0, posY=0;
    
    if(n==1)
      return new LineFunction(true);
    
    for(int i=0; i<n; i++){
      xSum+=((Point *)pointList[i])->x;
      ySum+=((Point *)pointList[i])->y;
      x2+=pow((double)(((Point *)pointList[i])->x), 2);
      xy+=((Point *)pointList[i])->x * ((Point *)pointList[i])->y;
    }
    
    double m=(double)(n*xy - xSum*ySum)/(double)(n*x2 - pow((double)xSum, 2) + 0.0000001);
    double b=(double)(ySum - m*xSum)/(double)n;
    return new LineFunction(m, b);        
  }
  
  
  void connectPoints(Point *p1, Point *p2, double *skel, int width, int height){
    double slope = (double)(p1->y - p2->y)/(double)(p1->x - p2->x);
    double b = p1->y - slope*p1->x;
    if(abs(p1->x - p2->x) > abs(p1->y - p2->y)){
      int begin = min(p1->x, p2->x), end = max(p1->x, p2->x);
      for(int x=begin; x <= end; x++){
	int y = ((int)(slope*x+b));
	if(y<0 || y>=height || x<0 || x>=width)
	  continue;
	skel[y*width + x] = -2;
      }
    }
    else{
      int begin = min(p1->y,p2->y), end = max(p1->y,p2->y);
      for(int y=begin; y <= end; y++){
	int x = (int)((y-b)/slope);
	if(x<0 || x>=width || y<0 || y>=height)
	  continue;
	skel[y*width + x] = -2;
      }
    }
  }
  
  bool elementOf(vector<Point *> v, Point *p){
    for(int i=0; i<v.size(); i++){
      Point *tmp = (Point *)v[i];
      if(p->x==tmp->x && p->y==tmp->y)
        return true;
    }
    return false;
  }

  double *thresh(double *image, double thresh, int width, int height){
    double *thresholdedImage = new double[width*height];
    for(int i=0; i < width; i++)
      for(int j=0; j < height; j++)
        thresholdedImage[j*width+i] = (image[j*width+i] > thresh) ? 1 : 0;
    
    return thresholdedImage;
  }
  
  
 public:
  Morphology(){
    pi = 3.141592654;
  }
  
  double *doubleThreshold(double *image, double lowThresh, double highThresh, int width, int height){
//    double *lowThresholded = new double[width*height];
    double *lowThresholded = thresh(image, lowThresh, width, height);
//    double *highThresholded = new double[width*height];
    double *highThresholded = thresh(image, highThresh, width, height);
    double *thresholdedImage = new double[width*height];
    queue<Point> nextPoints;

    bool *visited = new bool[width*height];
    for(int i=0; i < width*height; i++) {
      thresholdedImage[i] = 0;
      visited[i] = false;
    }
    
    // This alorithm seems to look for a lit pixel on the highthreshold image 
    // and then see if any of its neighbors are in the low threshold
    // I think the assumption is that the high-threshold pixels will always be
    // on an edge, and its neighbors on the low-threshold image will be too
    for(int i=0; i < width; i++) {
      for(int j=0; j < height; j++) {
        if(visited[j*width+i])
          continue;
        if(highThresholded[j*width+i] != 0){
          //new Point(i, j)
//          Point p = {i, j};
          nextPoints.push(Point{i, j});
          thresholdedImage[j*width+i] = 1;
          do {
//            Point *tmp = nextPoints.front();
//            int i1 = tmp->x;
//            int j1 = tmp->y;
//            nextPoints.pop();
            Point tmp = nextPoints.front();
            int i1 = tmp.x;
            int j1 = tmp.y;
            nextPoints.pop();
            visited[j1*width+i1] = true;
            // This checks to see if the neighbors are in the low threshold
            // First seeing if the neighbor is on the image
            for(int m=-1; m < 2; m++) {
              for(int n=-1; n < 2; n++) {
                if(i1+m < 0 || i1+m >= width || j1+n < 0 || j1+n >= height)
                    continue;
                if(!(n==0 && m==0) && !visited[((j1+n) * width) + i1 + m]){
                  if(lowThresholded[((j1+n) * width) + (i1+m)] != 0){
                    thresholdedImage[((j1+n) * width) + (i1+m)] = 1;
                    visited[((j1+n) * width) + (i1+m)] = true;
                    Point p = {i1 + m, j1 + n};
                    nextPoints.push(p);
                  }
                }

              }
            }

          }
          while(!nextPoints.empty());
        }
      }
    }
    
    delete[] visited;
    delete[] lowThresholded;
    delete[] highThresholded;
    return thresholdedImage;
  }


  
  double *detangle(double *image, int width, int height){
    printf("Detangling Image...\n");
    
    double threshold = (double)(width*height)/1200;
    tripleDetector *td = new tripleDetector();
    vector<Point *> triples = td->detectTriples(image, width, height);
    double *detangledImage { new double[width * height]{} };
    bool *visited{ new bool[width * height]{} };
    vector<Point *> region;
    double *intermediateImage{ new double[width*height]{} };
    //double *intermediateImage = new double[width*height];
    /*
//    double *detangledImage = new double[width * height];
//    bool *visited = new bool[width * height];
    for (int k = 0; k < width * height; k++){
      detangledImage[k] = 0;
      visited[k] = false;
    }
    */
    queue<Point *> nextPoints;

//    printf("TestPoint0\n");
    
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
//        printf("%f\n", image[j * width + i]);
        // if a pixel is 0 and hasn't been visited
        if (image[j * width + i] == 0.0 && !visited[j * width + i]) {
//          printf("in");
          int i1, j1;
          nextPoints.push(new Point(i, j));
          do{
            Point *tmp = (Point *)nextPoints.front();
            nextPoints.pop();
            i1 = tmp->x;
            j1 = tmp->y;
            if(region.size() <= threshold)
              region.push_back(new Point(i1, j1));

            visited[j1 * width + i1] = true;
            for (int x = -1; x <= 1; x++) {
              for (int y = -1; y <= 1; y++){
                if (x==y || x==-y || y==-x)
                  continue;
                int nextI = i1 + x;
                int nextJ = j1 + y;
                if (nextJ >= 0 && nextJ < height && nextI >= 0 && nextI < width){
                  if (!visited[nextJ * width + nextI]){
                    visited[nextJ * width + nextI] = true;
                    if (image[nextJ * width + nextI] == 0)
                      nextPoints.push(new Point(nextI, nextJ));
                  }
                }
              }
            }
          } while (!nextPoints.empty());
          
          if(region.size() < threshold){
            vector<Point *> regionBoundary;
            for(int k=0; k<region.size(); k++){
              Point *tmp = (Point *)region[k];
              i1 = tmp->x;
              j1 = tmp->y;
              for(int m=-1; m<=1; m++){
                for(int n=-1; n<=1; n++){
                  if((m==0 && n==0) || (i1+m<0 || i1+m >= width || j1+n<0 || j1+n>=height))
                    continue;
                  int i2 = i1+m;
                  int j2 = j1+n;
                  if(image[j2*width+i2] != 0)
                    regionBoundary.push_back(new Point(i2, j2));
                }
              }
            }
            
            vector<Point *> endPoints;
            for(int g=0; g<regionBoundary.size(); g++){
              Point *tmp = (Point *)regionBoundary[g];
              i1 = tmp->x;
              j1 = tmp->y;
              
              detangledImage[j1*width+i1]=-1;
              
              bool endLoops=false;
              for(int m=-1; m<=1; m++){
                for(int n=-1; n<=1; n++){
                  int i2=i1+m;
                  int j2 = j1+n;
                  if((m==0 && n==0) || (i2<0 || i2>=width || j2<0 || j2>=height) || elementOf(regionBoundary, new Point(i2, j2)) || elementOf(endPoints, new Point(i2, j2)))
                    continue;
                  if(image[j2*width+i2] != 0){
                    endPoints.push_back(new Point(i2, j2));
                    endLoops=true;
                    break;
                  }
                }
                if(endLoops)
                  break;
              }
            }
          
//l            double *intermediateImage = new double[width*height];
            for(int i=0; i<width*height; i++)
              intermediateImage[i] = (detangledImage[i]==-1)?0:image[i];
            
            if(endPoints.size() == 2)
              connectPoints(endPoints[0], endPoints[1], detangledImage, width, height);
            else if(endPoints.size() == 3){ //compute where the endpoints, when extended from the skeleton, will intersect each other, connect all endpoints to the intersection point that minimizes the angle between the point that is not part of the intersection and the intersection point, and the line emanating from the extra point
              int cx=0, cy=0;
              
              LineFunction *lf1 = functionOfEndPoint(intermediateImage, (Point *)endPoints[0], width, height, triples);
              LineFunction *lf2 = functionOfEndPoint(intermediateImage, (Point *)endPoints[1], width, height, triples);
              LineFunction *lf3 = functionOfEndPoint(intermediateImage, (Point *)endPoints[2], width, height, triples);
              
              if(!(lf1->bad || lf2->bad || lf3-> bad)){
            
                double m1=lf1->m, m2=lf2->m, m3=lf3->m;
                double b1=lf1->b, b2=lf2->b, b3=lf3->b;
                
                int x1=endPoints[0]->x, x2=endPoints[1]->x, x3=endPoints[2]->x;
                int y1=endPoints[0]->y, y2=endPoints[1]->y, y3=endPoints[2]->y;
                
                Point *intersection1 = (m1==m2) ? new Point((int)((x1+x2)/2), (int)((y1+y2)/2)) : new Point((int)((b2-b1)/(m1-m2)), (int)(m1*(b2-b1)/(m1-m2)+b1));
                Point *intersection2 = (m1==m3) ? new Point((int)((x1+x3)/2), (int)((y1+y3)/2)) : new Point((int)((b3-b1)/(m1-m3)), (int)(m1*(b3-b1)/(m1-m3)+b1));
                Point *intersection3 = (m3==m2) ? new Point((int)((x3+x2)/2), (int)((y3+y2)/2)) : new Point((int)((b2-b3)/(m3-m2)), (int)(m3*(b2-b3)/(m3-m2)+b3));
                
                int ix1=intersection1->x, ix2=intersection2->x, ix3=intersection3->x;
                int iy1=intersection1->y, iy2=intersection2->y, iy3=intersection3->y;
                
                double mn1 = (double)(iy1 - y3) / (double)(ix1 - x3);
                double mn2 = (double)(iy2 - y2) / (double)(ix2 - x2);
                double mn3 = (double)(iy3 - y1) / (double)(ix3 - x1);
                
                double angle1=myAtanDiff(m3, mn1), angle2=myAtanDiff(m2, mn2), angle3=myAtanDiff(m1, mn3);
                double minVal=min(angle1, min(angle2, angle3));
                
                Point *connectionPt;
                if(angle1 == minVal)
                  connectionPt = new Point(intersection1);
                else if(angle2 == minVal)
                  connectionPt = new Point(intersection2);
                else
                  connectionPt = new Point(intersection3);
                
                
                double maxDistance = 1.2*max(((Point *)endPoints[0])->distanceFrom(endPoints[1]), max(((Point *)endPoints[0])->distanceFrom(endPoints[2]), ((Point *)endPoints[1])->distanceFrom(endPoints[2])));
                if(!(connectionPt->distanceFrom(endPoints[0]) > maxDistance || connectionPt->distanceFrom(endPoints[1]) > maxDistance ||connectionPt->distanceFrom(endPoints[2]) > maxDistance)){
                  connectPoints(endPoints[0], connectionPt, detangledImage, width, height);
                  connectPoints(endPoints[1], connectionPt, detangledImage, width, height);
                  connectPoints(endPoints[2], connectionPt, detangledImage, width, height);
                }
              }
              
            }
            
            
          }
          region.clear();
        }
//        printf("%d, %d\n", i, j);
	
      }
//      printf(".");
    }

    delete[] intermediateImage;
    delete[] visited;

    for(int i=0; i<width*height; i++)
      detangledImage[i] = (detangledImage[i] == -1) ? 0 : ((detangledImage[i]==-2) ? 1 : image[i]);
    
    return detangledImage;
  }





  

  double *dilate(double *image, int width, int height){
    double *dilatedImage = new double[width * height];
    for(int i=0; i<width*height; i++)
      image[i] = (image[i]!=0) ? 1 : 0;
    for (int i = 1; i < width - 1; i++)
      for (int j = 1; j < height - 1; j++)
	if (image[j * width + i] != 0)
	  for (int b = -1; b <= 1; b++)
	    for (int n = -1; n <= 1; n++)
	      dilatedImage[(j + n) * width + (i + b)] = 1;
    return dilatedImage;
  }


  double *binaryDenoise(double *image, int width, int height, int threshold, int neighThresh){
    printf("Binary Denoising Image...\n");
    double *denoisedImage = new double[width * height];
    bool *visited = new bool[width * height];
    vector<Point> region;

    for (int k = 0; k < width * height; k++){
      denoisedImage[k] = 0;
      visited[k] = false;
    }

    queue<Point> nextPoints;

    // Loop through every pixel, and if its 0 in the denoised image an not zero in the normal image
    //
    for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++){
        if (denoisedImage[j * width + i] == 0 && image[j * width + i] != 0) {
          int i1, j1;
          Point p = {i, j};
          nextPoints.push(p);
          do {
            Point tmp = (Point) nextPoints.front();
            nextPoints.pop();
            i1 = tmp.x;
            j1 = tmp.y;
//            Point *tmp = (Point *)nextPoints.front();
//            nextPoints.pop();
//            i1 = tmp->x;
//            j1 = tmp->y;
            if (region.size() < threshold) {
              Point tmpr = {i1, j1};
              region.push_back(tmpr);
            }
            else if (region.size() == threshold) {
              for (int h = 0; h < region.size(); h++){
                Point p = region.at(h);
                denoisedImage[p.y * width + p.x] = 1;
              }
              denoisedImage[j1 * width + i1] = 1;
            } 
            else
              denoisedImage[j1 * width + i1] = 1;
            
            visited[j1 * width + i1] = true;

            for (int x = -1 * neighThresh; x <= neighThresh; x++)
              for (int y = -1 * neighThresh; y <= neighThresh; y++){
                if (x == 0 && y == 0)
                  continue;
                int nextI = i1 + x;
                int nextJ = j1 + y;
                // Checking if the neighbor is on a border of an image
                if (nextJ >= 0 && nextJ < height && nextI >= 0 && nextI < width)
                  if (!visited[nextJ * width + nextI]) {
                    visited[nextJ * width + nextI] = true;
                    if (image[nextJ * width + nextI] != 0 && denoisedImage[nextJ * width + nextI] == 0)
                      Point p = {nextI, nextJ};
                      nextPoints.push(p);
                  }
              }
          }
          while (!nextPoints.empty());
          region.clear();
        }
      }
    delete[] visited;
    return denoisedImage;
  }


  


  double *thin(double *boundaryMap, int width, int height){

    // Louisa Lam, Seong-Whan Lee, and Ching Y. Wuen, "Thinning Methodologies-A
    // Comprehensive Survey," IEEE TrPAMI, vol. 14, no. 9, pp. 869-885, 1992.  The
    // algorithm is described at the bottom of the first column and the top of the
    // second column on page 879.


    bool G1[] = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    1,    1,    1,    1,    0,    1,    1,
		 1,    1,    1,    1,    0,    0,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    0,    1,    1,    1,    0,    1,    1,    0,    0,    1,    1,
		 0,    0,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    0,    0,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    1,    1,    0,    0,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 0,    0,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    1,    1,
		 1,    0,    1,    1,    1,    1,    0,    0,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    0,    0,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    0,    1,    1,    1,    0,    1,    1,    1,    1,    0,    0,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    0,    1,    1,    1,    0,    1,    1,
		 0,    0,    1,    1,    0,    0,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    0,    0,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    0,    1,    1,    1,    0,    1,    1,    1,    1,    0,    0,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    0,    0,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    0,    1,    1,    1,    0,    1,    1,
		 1,    1,    0,    0,    1,    1,    0,    0};
    
    bool G2[] = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    1,
		 0,    0,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    1,    1,    1,    0,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,
		 1,    1,    1,    1,    0,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    0,    0,    1,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    0,    1,    0,    1,    0,    1,    1,    0,    0,
		 1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    1,    0,    1,    0,    1,    1,    1,    1,
		 1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    0,
		 1,    0,    1,    0,    1,    1,    0,    0,    1,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    0,    0,
		 1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 1,    0,    1,    0,    1,    1,    1,    1,    1,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    0,    1,    0,    1,    0,
		 1,    1,    0,    0,    1,    0,    0,    0};
    
    bool G3[] = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,
		 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    1,    1,    1,    1,    1,    1,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0};
    
    bool G4[] = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 1,    1,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    1,
		 1,    1,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    0,    1,    1,    1,    0,    0,    1,    1,    0,    1,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 1,    1,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    1,
		 1,    1,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    0,    1,    1,    1,    0,    0,    1,    1,    0,    1,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 1,    1,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    1,
		 1,    1,    0,    0,    1,    1,    0,    1,    1,    1,    0,    0,
		 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		 0,    0,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,
		 1,    1,    0,    0,    1,    1,    0,    0};
    
    int size = width * height;
    bool *a = new bool[size];
    bool *b = new bool[size];
    bool *c = new bool[size];
    bool *d = new bool[size];

  
    for (int i = 0; i < size; i++)
      a[i] = (bool)((boundaryMap[i] != 0) ? 1 : 0);
      
    for(int i=0; i<size; i++){
      b[i] = a[i];
      c[i] = a[i];
    }
    
    int *lookup = new int[size];
    bool *G1lookup = new bool[size];
    bool *G2lookup = new bool[size];
    bool *G3lookup = new bool[size];
    bool *G4lookup = new bool[size];
    while (true){
      // Make a lookup table that will produce
      // a lookup table indices. This is avoid
      // doing more work in calling applylut
      // multiple times with the same input than
      // we really need to.
      // >> lutlut = 1:512;
      int *lutlut = new int[512];
      for (int i = 0; i < 512; i++)
        lutlut[i] = i;

      // Apply the lutlut LUT to a, yielding lookup
      // >> lookup = applylut(a, lutlut);
//      int *lookup = new int[size];
      for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
          lookup[x + y * width] = lutlut[Nhood3Offset(c, width, height, x, y)];

      // First subiteration
      // >> d(:) = G1(lookup) & G2(lookup) & G3(lookup);
//      bool *G1lookup = new bool[size];
//      bool *G2lookup = new bool[size];
//      bool *G3lookup = new bool[size];
//      bool *G4lookup = new bool[size];
      for (int i = 0; i < size; i++){
        G1lookup[i] = G1[lookup[i]];
        G2lookup[i] = G2[lookup[i]];
        G3lookup[i] = G3[lookup[i]];
      }
      for (int i = 0; i < size; i++)
        d[i] = (bool)(G1lookup[i] & G2lookup[i] & G3lookup[i]);
      // >> c = a & ~d;
      for (int i = 0; i < size; i++)
        c[i] = (bool)(c[i] & ~d[i]);

      // Second subiteration
      // >> lookup = applylut(c, lutlut);
//      G1lookup = new bool[size];
//      G2lookup = new bool[size];
//      G3lookup = new bool[size];
//      G4lookup = new bool[size];
//      lookup = new int[size];

      for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
          lookup[x + y * width] = lutlut[Nhood3Offset(c, width, height, x, y)];
      // >> d(:) = G1(lookup) & G2(lookup) & G4(lookup);
      for (int i = 0; i < size; i++){
        G1lookup[i] = G1[lookup[i]];
        G2lookup[i] = G2[lookup[i]];
        G4lookup[i] = G4[lookup[i]];
      }
      for (int i = 0; i < size; i++)
        d[i] = (bool)(G1lookup[i] & G2lookup[i] & G4lookup[i]);
      // >> c = c & ~d;
      for (int i = 0; i < size; i++)
        c[i] = (bool)(c[i] & ~d[i]);

      bool differs = false;
      for (int i = 0; i < size; i++)
        if (c[i] != b[i])
          { differs = true; break; }
      if (!differs) break;

      // Remember old result for next cycle
      for(int i=0; i<size; i++)
        b[i] = c[i];

    }

    delete[] lookup;
    delete[] G1lookup;
    delete[] G2lookup;
    delete[] G3lookup;
    delete[] G4lookup;

    double *resultMapBitmap = new double[size];
    for(int i=0; i<size; i++)
      resultMapBitmap[i] = (double)c[i];

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    return resultMapBitmap;
  }

  int Nhood3Offset(bool *logicalData, int width, int height, int x, int y){
    int weights3[3][3];
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        weights3[i][j] = pow(2.0, i*3+j);

    int minX, maxX, minY, maxY;
    int result = 0;

    // Figure out the neighborhood extent that does not go past image boundaries
    minX = ((x == 0) ? (1) : (0));
    maxX = ((x == (width - 1)) ? (1) : (2));
    minY = ((y == 0) ? (1) : (0));
    maxY = ((y == (height - 1)) ? (1) : (2));

    for (int xx = minX; xx <= maxX; xx++)
      for (int yy = minY; yy <= maxY; yy++)
	if (logicalData[(x + xx - 1) + (y + yy - 1) * width] != 0)
	  result += weights3[xx][yy];

    return result;
  }


/*
%
% Function THIN
%
function [c,lut] = thin(a)
  
% Louisa Lam, Seong-Whan Lee, and Ching Y. Wuen, "Thinning Methodologies-A
% Comprehensive Survey," IEEE TrPAMI, vol. 14, no. 9, pp. 869-885, 1992.  The
% algorithm is described at the bottom of the first column and the top of the
% second column on page 879.

*/




  double *prune(double *skel, int threshold, int width, int height){
    printf("Pruning Image...\n");
  
//    printf("Testpoint0\n");
    double *clippedImage = new double[width*height];
    
    for(int i=0; i<width*height; i++)
      clippedImage[i]=skel[i];
  
    cornerDetector *cd = new cornerDetector();
    vector<Point *> corners = cd->detectCorners(skel, width, height);
    endDetector *ed = new endDetector();
    vector<Point *> ends = ed->detectEnds(skel, width, height);
    tripleDetector *td = new tripleDetector();
    vector<Point *> trips = td->detectTriples(skel, width, height);
    
    bool *interestPts = new bool[width*height];
    for(int i=0; i<width*height; i++)
      interestPts[i]=false;
    for(int i=0; i<trips.size(); i++){
      Point *tmp = trips[i];
      interestPts[tmp->y*width+tmp->x]=true;
    }
  
//    printf("Testpoint1\n");
    for(int i=0; i<ends.size(); i++){
      Point *tmp = ends[i];
      bool done=false;
      bool foundNeighbor=true;
      int count=0;
      vector<Point *> segment;
      segment.push_back(new Point(tmp->x, tmp->y));
      while(!done && foundNeighbor) {
        foundNeighbor=false;
        for(int x=tmp->x-1; x<=tmp->x+1; x++) {
          for(int y=tmp->y-1; y<=tmp->y+1; y++) {
            if (x<0 || x>width-1 || y<0 || y>height-1 || (y==tmp->y && x==tmp->x))
              continue;
            if (skel[y*width+x]!=0  && !elementOf(segment, new Point(x, y))){
              segment.push_back(new Point(x, y));
              tmp->x=x;
              tmp->y=y;
              foundNeighbor=true;
              goto hasNeighbor;
            }
            if (interestPts[y*width+x])
              done=true;
          }
        }
        hasNeighbor:; // the semi-colon is an empty statement
      }

      if(segment.size() < threshold)
        for(int k=0; k<segment.size(); k++) {
          Point *p = (Point *)segment[k];
          clippedImage[p->y*width+p->x]=0;
        }
    }
    
    delete[] interestPts;
    return clippedImage;
    
  }



};
