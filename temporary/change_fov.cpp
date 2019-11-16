#include <iostream>
#include "fileio.h"


int main()
{
  int w = 512;
  int h = 256;
  int d = 180;
  std::vector<float> f(w * h * d, 0.);
  std::vector<float> g(f.size(), 0.);
  std::vector<int> fov(w * h, 0.);
  readRawFile("./img/primary_detector_float_512-256-180.raw", f);
  readRawFile("./img/fov_11pinhole_5mm_int_512-256.raw", fov);
  for(int k = 0; k < d; k++)
  {
    for(int i = 0; i < h; i++)
    {
      for(int j = 0; j < w; j++)
      {
        bool remove = (fov[i * w + j] == 1 || fov[i * w + j] == 2 || fov[i * w + j] == 4 || fov[i * w + j] == 6 || fov[i * w + j] == 8 || fov[i * w + j] == 9 );
        int index = k * w * h + i * w + j;
        g[index] = remove ? 0. : f[index];
      }
    }
  }
  writeRawFile("remove_primary_detector_float_512-256-180.raw", g);

}
