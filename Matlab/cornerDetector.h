#ifndef RAY
#define RAY
#include "Ray.h" 
#endif 
#ifndef POINT
#define POINT
#include "Point.h" 
#endif 
#include <math.h>
#include <vector>
#include <queue>
#define numPointsForSlope 15
#define cornerThresh 2.4 //smaller values mean fewer corners
#define clusterThresh 4
using namespace std;



class PointDiff {
 public:
  Point *p;
  Point *ngh1;
  Point *ngh2;
  double angle;
 
  PointDiff(Point *point, double a, Point *n1, Point *n2){
    p = new Point(point);
    angle = a;
    ngh1 = new Point(n1);
    ngh2 = new Point(n2);
  }

  PointDiff(Point *point, double a){
    p = new Point(point);
    angle = a;
    ngh1 = NULL;
    ngh2 = NULL;
  }

  PointDiff(PointDiff *pd){
    p = new Point(pd->p);
    angle = pd->angle;
    ngh1 = new Point(pd->ngh1);
    ngh2 = new Point(pd->ngh2);
  }
};




class cornerDetector{
 private:

  double pi;

  bool elementOf(vector<Point *> v, Point *p){
    for(int i=0; i<v.size(); i++){
      Point *tmp = (Point *)v[i];
      if(p->x==tmp->x && p->y==tmp->y)
	return true;
    }
    return false;
  }
  

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



  Ray *ray(Point *p, double *skel, int i, int j, int width, int height, vector<Point *> trips){
  vector<Point *> pointList;
  Point *tripPoint=NULL;
  Point *point = new Point(p);
  pointList.push_back(new Point(point));
  while(pointList.size()<numPointsForSlope && point!=NULL) {
    int x=point->x;
    int y=point->y;
    point=NULL;
    bool breakLoops=false;
    for(int m=x-1; m<=x+1; m++) {
      for(int n=y-1; n<=y+1; n++) {
        if((m==i && n==j) || (m==x && n==y) || m<0 || m>=width || n<0 || n>=height)
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
  int xDirection, yDirection;
  
  if(n==1){
    pointList.clear();
    pointList.push_back(new Point(i, j));
    pointList.push_back(new Point(p));
    if(tripPoint!=NULL)
      pointList.push_back(new Point(tripPoint));
    n = pointList.size();
  }

  for(int i=0; i<n; i++){
    xSum+=((Point *)pointList[i])->x;
    ySum+=((Point *)pointList[i])->y;
    x2+=pow((double)(((Point *)pointList[i])->x), 2);
    xy+=((Point *)pointList[i])->x * ((Point *)pointList[i])->y;
    if(i!=n-1){
      if( ( (int)((Point *)pointList[i])->x - (int)((Point *)pointList[i+1])->x) < 0)
        posX++;
      else if( ( (int)((Point *)pointList[i])->x - (int)((Point *)pointList[i+1])->x) > 0)
        negX++;
      if( ( (int)((Point *)pointList[i])->y - (int)((Point *)pointList[i+1])->y) < 0)
        posY++;
      else if( ( (int)((Point *)pointList[i])->y - (int)((Point *)pointList[i+1])->y) > 0)
        negY++;
    }
  }
  xDirection = (negX>posX) ? -1 : 1;
  yDirection = (negY>posY) ? -1 : 1;

  if(n*x2==pow((double)xSum, 2))
    return new Ray(1000, xDirection, yDirection);
  
  return new Ray((double)(n*xy - xSum*ySum)/(double)(n*x2 - pow((double)xSum, 2) + 0.000001), xDirection, yDirection);
}


  double absoluteAngle(Ray *ray){
    int q = ray->quad;
    double angle = ray->refAngle();
    if(q==1)
      return angle;
    if(q==2)
      return pi-angle;
    if(q==3)
      return pi+angle;
    return 2*pi-angle;
  }

  
  double rayDiff(Ray *r1, Ray *r2){
    double absAngle1 = absoluteAngle(r1), absAngle2 = absoluteAngle(r2);
    return min(abs(absAngle1 - absAngle2), abs(2*pi - abs(absAngle1 - absAngle2)));
  }


 public:

  cornerDetector(){
    pi = 3.1415926535897932384626433832795;
  }

  vector<Point *> detectCorners(double *skeleton, int width, int height) {
    vector<Point *> corners;
    vector<PointDiff *> possibleCorners;
    vector<Point *> trips;
    for(int i=0; i<width; i++) {
      for(int j=0; j<height; j++) {
        if(skeleton[j*width+i]!=0) {
          vector<Point *> neighbors;
          for(int m=-1; m<=1; m++) {
            for(int n=-1; n<=1; n++) {
              int x = i+m, y = j+n;
              if( (x==i && y==j) || x<0 || x>=width || y<0 || y>=height)
                continue;
              if(skeleton[y*width+x] != 0)
                neighbors.push_back(new Point(x, y));
            }
          }
          int nghSize = neighbors.size();
          if(nghSize >= 2) {
            Point *ngh1 = new Point((Point *)neighbors[0]);
            Point *ngh2 = NULL;
            int numClusters = nghSize;
            for(int k=0; k<nghSize; k++){
              Point *pt1 = new Point((Point *)neighbors[k]);
              if(pt1->distanceFrom(ngh1) > 1)
                ngh2 = new Point(pt1);
              for(int h=k; h<nghSize; h++){
                Point *pt2 = new Point((Point *)neighbors[h]);
                if(pt1->distanceFrom(pt2)==1)
                  numClusters--;
              }
            }
            if(numClusters==2 && ngh2!=NULL)
              possibleCorners.push_back(new PointDiff(new Point(i, j), 0, new Point(ngh1), new Point(ngh2)));
            else if(numClusters>2)
              trips.push_back(new Point(i, j));
          }
          
        }
      }
    }
      
    int n = possibleCorners.size();
    vector<PointDiff *> refinedCornerList;
        
    for(int i=0; i<n; i++){
      Point *p1 = (Point *)((PointDiff *)possibleCorners[i])->p;
      Point *ngh1 = (Point *)((PointDiff *)possibleCorners[i])->ngh1;
      Point *ngh2 = (Point *)((PointDiff *)possibleCorners[i])->ngh2;
      int x = p1->x;
      int y = p1->y;
      if(elementOfTrip(trips, new Point(x,y), 2))
        continue;
      Ray *ray1 = ray(ngh1, skeleton, x, y, width, height, trips);
      Ray *ray2 = ray(ngh2, skeleton, x, y, width, height, trips);
      double angle = rayDiff(ray1, ray2);
      if(angle >= cornerThresh || angle <= pi-cornerThresh)
        continue;
      refinedCornerList.push_back(new PointDiff(new Point(x, y), angle));
      /*************************************/
      //corners.push_back(new Point(x, y));
      /*************************************/
    }
    printf("tp0.1");
    
    double *cornerAngles = new double[width*height];
    bool *used = new bool[width*height];
    for(int i=0; i<width*height; i++){
      cornerAngles[i]=0;
      used[i]=false;
    }
    for(int i=0; i<refinedCornerList.size(); i++){
      PointDiff *tmp = (PointDiff *)refinedCornerList[i];
      Point *pt = tmp->p;
      cornerAngles[pt->y*width+pt->x] = (double)tmp->angle;
    }
    for(int i=0; i<width; i++) {
      for(int j=0; j<height; j++) {
        if(cornerAngles[j*width+i]!=0 && !used[j*width+i]){
          Point *nextPoint = new Point(i, j);
          vector<PointDiff *> pointCluster;
          while(nextPoint != NULL){
            int x = nextPoint->x;
            int y = nextPoint->y;
            delete nextPoint;
            nextPoint = NULL;
            used[y*width+x]=true;
            bool breakLoops=false;
            for(int m=x-clusterThresh; m<=x+clusterThresh; m++){
              for(int n=y-clusterThresh; n<=y+clusterThresh; n++){
                if((m == x && n == y) || m<0 || n<0 || m>=width || n>=height)
                  continue;
                if(cornerAngles[n*width+m]!=0 && !used[n*width+m]){
                  pointCluster.push_back(new PointDiff(new Point(m, n), cornerAngles[n*width+m]));
                  nextPoint = new Point(m, n);
                  breakLoops=true;
                }
                used[n*width+m]=true;
                if(breakLoops)
                  break;
              }
              if(breakLoops)
                break;
            }
          }
          double minAngle=10000000;
          Point *minPoint;
          for(int k=0; k<pointCluster.size(); k++){
            PointDiff *tmpDiff = pointCluster[k];
            if(tmpDiff->angle < minAngle){
              minAngle = tmpDiff->angle;
              minPoint = new Point(tmpDiff->p);
            }
          }
          corners.push_back(new Point(minPoint));
        }
      }
    }

    return corners;
  }
  
  
};
