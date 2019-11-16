#include <iostream>
#include <vector>
#include "fileio.h"

int main()
{
  std::vector<float> f(512 * 256);
  std::vector<int> g(f.size(), -1);
  readRawFile("fov_float_512-256.raw", f);
  // readRawFile("primary_single_pinhole_float_512-256-180.raw", f);

  // for(int i = 0; i < 512; i++)
  // {
  //   for(int j = 0; j < 256; j++)
  //   {
  //     g[i * 512 + j] = (int)(f[i * 512 + j]);
  //   }
  // }

  for(int i = 0; i < f.size(); i++) g[i] = (int)(f[i]);

  writeRawFile("fov_int_512-256.raw", g);
  // writeRawFile("proj_n_float_512-180.raw", g);

}
