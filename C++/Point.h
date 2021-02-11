//Code modified 6/30/20 by Ryan Cecil
//Operators did not need to be defined by scope since
//they were located in the class.

#include <math.h>
#include <stdio.h>

//Construct clas Point
class Point{
 public:
	int x;
	int y;
	bool endpoint;
	Point(int x1, int y1){
		x=x1;
		y=y1;
		endpoint=false;
	}
	Point(Point *p){
		x=p->x;
		y=p->y;
		endpoint=false;
	}
	Point(int x1, int y1, bool endp){
		x=x1;
		y=y1;
		endpoint=endp;
	}
	Point(Point *p, bool endp){
		x=p->x;
		y=p->y;
		endpoint=endp;
	}
	Point(Point p, bool endp){
		x=p.x;
		y=p.y;
		endpoint=endp;
	}
	Point(){
		endpoint=false;
	}
	double distanceFrom(Point *p){
		return sqrt(pow(x-(p->x), 2.0) + pow(y-(p->y), 2.0));
	}
	bool equals(Point *p){
		return (x==p->x && y==p->y) ? true : false;
	}
	char *toString(){
	  char *pt = new char[100];
	  sprintf(pt, "(%d, %d)", x, y);
	  return pt;
    }
      
    //Simply removed the Point Point:: part of the code.
    Point operator+= (Point *b) { 
    x+=b->x;
    y+=b->y;
    return this;  
    }
    Point operator+ (Point *b) { 
        return new Point(x+ b->x, y+b->y);  
    }
    Point operator/ (Point *b) { 
        return new Point((int)(x/b->x), int(y/b->y));  
    }
    Point operator/ (int b) { 
        return new Point((int)(x/b), int(y/b));  
    }
};

