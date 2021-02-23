#ifndef POINT
#define POINT
#include "Point.h" 
#endif 

#include <stdlib.h>
#include <vector>
using namespace std;


class tripleDetector{
 public:
  vector<Point *>detectTriples(double *skeleton, int width, int height){
    vector<Point *> triples;
    for(int i=0; i<width; i++)
      for(int j=0; j<height; j++){
	if(skeleton[j*width+i] != 0){
	  vector <Point *> neighbors;
	  for(int m=-1; m<=1; m++)
	    for(int n=-1; n<=1; n++){
	      int x = i+m, y = j+n;
	      if( (x==i && y==j) || x<0 || x>=width || y<0 || y>=height)
		continue;
	      if(skeleton[y*width+x] != 0)
		neighbors.push_back(new Point(x, y));
	    }
	  int nghSize = neighbors.size();
	  if(nghSize>2){
	    int numClusters=nghSize;
	    for(int k=0; k<nghSize; k++){
	      Point *pt1 = (Point *)neighbors[k];
	      for(int h=k; h<nghSize; h++){
		Point *pt2 = (Point *)neighbors[h];
		if(pt1->distanceFrom(pt2)==1)
		  numClusters--;
	      }
	    }
	    if(numClusters>2)
	      triples.push_back(new Point(i,j));
	  }
	  
	}
      }
    
    return triples;
  }
};
