#include <iostream>
#include <vector>
#include "fileio.h"


int main()
{
  std::vector<float> f(512*256, 0.);
  std::vector<float> g(f.size(), 0.);
  readRawFile("efficiency_map_float_512-256.raw", f);
  for(int i = 0; i < f.size(); i++)
  {
  	g[i] = f[i] > 0 ? 1. : -1.;

  }
  writeRawFile("fov_float_512-256.raw", g);

}
