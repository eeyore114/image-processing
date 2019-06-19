#include <iostream>
#include "fileio.h"


int main()
{
  std::vector<float> f(128*128*128, 0.);
  std::vector<float> g(f.size(), 0.);
  readRawFile("Shepp_float_128-128-128.raw", f);
  for(int i = 0; i < f.size(); i++) { g[i] = f[i] > 0 ? 0.15 : 0.; }
  writeRawFile("absorp_map_shepp_float_128-128-128.raw", g);

}
