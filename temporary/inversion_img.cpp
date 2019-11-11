#include <iostream>
#include "fileio.h"


int main()
{
  int w = 64;
  int h = 64;
  int d = 64;
  std::vector<float> f(w * h * d, 0.);
  std::vector<float> g(f.size(), 0.);
  std::vector<int> fov(w * h, 0.);
  readRawFile("./img/ML-EM100_float_64-64-64.raw", f);
  readRawFile("./img/fov_11pinhole_5mm_int_512-256.raw", fov);
  for(int k = 0; k < d; k++)
  {
    for(int i = 0; i < h; i++)
    {
      for(int j = 0; j < w; j++)
      {
        int index = k * w * h + i * w + j;
        g[index] = f[k * w * h + i * w + (w - j)];
      }
    }
  }
  writeRawFile("ML-EM100_invert_float_64-64-64.raw", g);

}
