#include <math.h>

class Ray {
 private:
  int quadrant(int x, int y){
    if(x == 1 && y == 1)
      return 1;
    if(x == -1 && y == 1)
      return 2;
    if(x == -1 && y == -1)
      return 3;
    if(x == 1 && y == -1)
      return 4;
    else
      return 0;
  }

 public:
  double slope;
  int xDirection;
  int yDirection;
  int quad;

  Ray(double s, int xDir, int yDir){
    slope = s;
    xDirection = xDir;
    yDirection = yDir;
    quad = quadrant(xDir, yDir);
  }

  Ray(double rise, double run){
    slope = rise/run;
    yDirection = (rise>0) ? 1 : -1;
    xDirection = (run>0) ? 1 : -1;
    quad = quadrant(xDirection, yDirection);
  }

  double refAngle(){
    return abs(atan(slope));
  }
};
