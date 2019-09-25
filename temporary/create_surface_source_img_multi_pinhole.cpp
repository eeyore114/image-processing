#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "fileio.h"

enum Coordinate { X, Y, Z, coord_num };

int main()
{
  int w = 2048;
  int h = 1024;
  float rotation_radius = 17;
  float pinhole_img_pixel_size = 0.02;
  std::vector<int> f(w * h, -1);
  std::vector<float> pinhole_center
	{
	    rotation_radius, - 11.,  4.,
	    rotation_radius, -  3., 4.5,
	    rotation_radius,    3., 4.5,
	    rotation_radius,   11.,  4.,
	    rotation_radius, - 6.5,  0.,
	    rotation_radius,    0.,  0.,
	    rotation_radius,   6.5,  0.,
	    rotation_radius, - 11., - 4.,
	    rotation_radius, -  3., - 4.5,
	    rotation_radius,    3., - 4.5,
	    rotation_radius,   11., - 4.
	};


  for(int i = 0; i < h; i++)
  {
    for(int j = 0; j < w; j++)
    {
      float x = (- (w - 1.0) / 2.0 + j) * pinhole_img_pixel_size;
    	float y =   ((h - 1.0) / 2.0 - i) * pinhole_img_pixel_size;
      for(int k = 0; k < 11; k++)
      {
        float scope = 0.8;
        float x1 = pinhole_center[coord_num * k + Y];
        bool in_scope_of_x = (x > x1 - scope && x < x1 + scope);
        float y1 = pinhole_center[coord_num * k + Z];
        bool in_scope_of_y = (y > y1 - scope && y < y1 + scope);
        if(in_scope_of_x && in_scope_of_y) f[i * w + j] = k;
      }
    }
  }
  writeRawFile("surfaceSource_float_2048-1024.raw", f);
}
