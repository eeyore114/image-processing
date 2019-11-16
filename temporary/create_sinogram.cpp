#include <iostream>
#include <vector>
#include "fileio.h"


int main()
{
  std::vector<float> f(512 * 256 * 180, 0.);
  std::vector<float> g(512 * 180, 0.);
  readRawFile("efficiency_correction_temporary_float_512-256-180.raw", f);
  // readRawFile("primary_single_pinhole_float_512-256-180.raw", f);

  for(int i = 0; i < 180; i++)
  {
    for(int j = 0; j < 512; j++)
    {
      g[i * 512 + j] = f[i * 512 * 256 + 127 * 512 + j];
    }
  }
  writeRawFile("proj2_float_512-180.raw", g);
  // writeRawFile("proj_n_float_512-180.raw", g);

}
