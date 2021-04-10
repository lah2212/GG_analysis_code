#ifndef POINT
#define POINT
#include "Point.h" 
#endif 

#ifndef RAY
#define RAY
#include "Ray.h" 
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

#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

class Interpolation{
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
  
  

  Ray *createRay(Point *p, double *skel, int i, int j, int width, int height, vector<Point *> trips){
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
    return min(fabs(absAngle1 - absAngle2), fabs(2*pi - fabs(absAngle1 - absAngle2)));
  }



  Ray *rayComposition(Ray *r1, Ray *r2){
    double absAngle, absAngle1 = absoluteAngle(r1), absAngle2 = absoluteAngle(r2);
    double theta = min(fabs(absAngle1 - absAngle2), 2*pi - fabs(absAngle1 - absAngle2));
    if(theta == fabs(absAngle1 - absAngle2))
      if(max(absAngle1, absAngle2) < pi)
        absAngle = max(absAngle1, absAngle2) + pi - theta/2;
      else
        absAngle = max(absAngle1, absAngle2) - pi - theta/2;
    else
      absAngle = max(absAngle1, absAngle2) - pi + theta/2;
    
    
    while(absAngle < 0 || absAngle > 2*pi){
      if(absAngle<0)
        absAngle+=2*pi;
      if(absAngle>2*pi)
        absAngle-=2*pi;
    }
    
    Point *quad;
    if(absAngle < pi/2)
      quad = new Point(1, 1);
    else if(absAngle < pi)
      quad = new Point(-1, 1);
    else if(absAngle < 3*pi/2)
      quad = new Point(-1, -1);
    else
      quad = new Point(1, -1);
    
    return new Ray(tan(absAngle), quad->x, quad->y);
    
  } 



  double f1(double s){
    return (s>1) ? 1/s : 1;
  }
 


  void drawConnection(double *skel, Point *p1, Point *p2, int width, int height){
    double slope = (double)(p1->y - p2->y)/(double)(p1->x - p2->x+0.00001);
    double b = p1->y - slope*p1->x;
    if(abs(p1->x - p2->x) > abs(p1->y - p2->y)){
      int begin = (p1->x < p2->x) ? p1->x : p2->x;
      int end = (p1->x > p2->x) ? p1->x : p2->x;
      for(int x=begin; x <= end; x++){
	int y = ((int)(slope*x+b));
	if(y<0 || y>=height || x<0 || x>=width)
	  continue;
	skel[y*width + x] = 1;
      }
    }
    else{
      int begin = (p1->y < p2->y) ? p1->y : p2->y;
      int end = (p1->y > p2->y) ? p1->y : p2->y;
      for(int y=begin; y <= end; y++){
	int x = (int)((y-b)/slope);
	if(x<0 || x>=width || y<0 || y>=height)
	  continue;
	skel[y*width + x] = 1;
      }
    }
  }
  
  
  

 public:
  Interpolation(){
    pi = 3.1415926535897932384626433832795;
  }
  
  double *interpolate(double *skeleton, int width, int height){
    
    double *interpolatedSkel = new double[width*height];
    for(int i=0;i<width*height;i++)
      interpolatedSkel[i]=0;
    
    
    double maxDist = sqrt((double)width*height) / 21;
    cornerDetector *cd = new cornerDetector();
    vector<Point *> corners = cd->detectCorners(skeleton, width, height);
    endDetector *ed = new endDetector();
    vector<Point *> ends = ed->detectEnds(skeleton, width, height);
    tripleDetector *td = new tripleDetector();
    vector<Point *> triples = td->detectTriples(skeleton, width, height);
    

    bool *connectedIndeces = new bool[ends.size() + corners.size()];
    for(int i=0; i<ends.size()+corners.size(); i++)
      connectedIndeces[i]=false;


    
    for(int i=0; i<ends.size(); i++){
      if(connectedIndeces[i])
      	continue;
      Point *sendPt = new Point(ends[i]);
      Ray *sendingRay = createRay(sendPt, skeleton, sendPt->x, sendPt->y, width, height, triples);

      vector<double> weights;
      vector<int> weightIndeces;

      for(int j=0; j < (ends.size() + corners.size()); j++){
	Point *receivePt = (j < ends.size()) ? new Point (ends[j]) : new Point(corners[j-ends.size()]);
	Ray *receivingRay;
	double distance = sendPt->distanceFrom(receivePt);
	if(distance > maxDist || i==j || connectedIndeces[j])
	  continue;


	if(j>=ends.size()){
	  int x = receivePt->x, y = receivePt->y;
	  Point *ngh1=NULL;
	  Point *ngh2=NULL;
	  for(int m=x-1; m<=x+1; m++){
	    for(int n=y-1; n<=y+1; n++){
	      if((m==x && n==y) || m<0 || n<0 || m>=width || n>=height){
            continue;
          }
          if(skeleton[n*width+m] != 0) {
              if(ngh1 == NULL){
                ngh1 = new Point(m, n);
              }
              else if(ngh1->distanceFrom(new Point(m, n)) > 1){
                ngh2 = new Point(m, n);
              }
          }
	    }
      }
	  if(ngh2==NULL){
	    continue;
	  }
	  Ray *ray1 = createRay(ngh1, skeleton, x, y, width, height, triples);
	  Ray *ray2 = createRay(ngh2, skeleton, x, y, width, height, triples);
	  receivingRay = rayComposition(ray1, ray2);
	}
	else
	  receivingRay = createRay(receivePt, skeleton, receivePt->x, receivePt->y, width, height, triples);


	//vector<Point *> connection = connectPoints(sendPt, receivePt);
	double averageGrad=0;
	//for(int k=0; k<connection.size(); k++) compute gradient under connecting line (make sure gradient is normalized between 0 and 1 - 0 near edge and 1 away)
	

	double distThresh = maxDist/5;
	double fd = f1(distance/distThresh);
	double wd = 0.2;
	double angleThresh = pi/12;
	Ray *connectingRay1 = new Ray((receivePt->y - sendPt->y), (receivePt->x - sendPt->x));
	Ray *connectingRay2 = new Ray((sendPt->y - receivePt->y), (sendPt->x - receivePt->x));
	double angle1 = rayDiff(sendingRay, connectingRay2);
	double angle2 = (j>=ends.size()) ? rayDiff(receivingRay, connectingRay2) : rayDiff(receivingRay, connectingRay1);
	if(angle1 > 2*pi/5 || angle2 > 2*pi/5)
	  continue;

	double fa1 = f1(angle1 / angleThresh);
	double fa2 = f1(angle2 / angleThresh);
	double gradThresh = 0.7;
	double wa = 0.25;
	double fg = f1(averageGrad / gradThresh);
	double wg = 0.3;
	
	double weight = wd*fd + wa*fa1 + wa*fa2 + wg*fg;
	weights.push_back(weight);
	weightIndeces.push_back(j);
      }


      double weightThreshold=0.2;
      double maxWeight=0;
      int maxIndex=-1;
      for(int b=0; b<weights.size(); b++)
	if(weights[b]>maxWeight){
	  maxWeight=weights[b];
	  maxIndex=weightIndeces[b];
	}
      if(maxWeight>weightThreshold){
	Point *point2 = (maxIndex < ends.size()) ? new Point(ends[maxIndex]) : new Point(corners[maxIndex - ends.size()]);
	drawConnection(interpolatedSkel, sendPt, point2, width, height);
	connectedIndeces[i]=true;
	connectedIndeces[maxIndex]=true;
      }
      
    }
    


    for(int i=0; i<width*height; i++)
      if(interpolatedSkel[i]==0)
	interpolatedSkel[i]=skeleton[i];

    return interpolatedSkel;
  }

};
      
