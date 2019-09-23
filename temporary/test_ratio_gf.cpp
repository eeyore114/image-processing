#include <iostream>
#include <vector>
#include "fileio.h"


int main()
{
  int w = 1024;
  int h = 512;
  int d = 180;
  std::vector<float> f(w * h * d, 0.);
  std::vector<float> g(f.size(), 0.);
  std::vector<float> ratio_gf(f.size(), 0.);
  readRawFile("./img/detec_map_float_1024-512-180.raw", f);
  readRawFile("./img/f_proj_float_1024-512-180.raw", g);
  // readRawFile("primary_single_pinhole_float_512-256-180.raw", f);

  // for(int i = 0; i < 180; i++)
  // {
  //   for(int j = 0; j < 512; j++)
  //   {
  //     g[i * 512 + j] = f[i * 512 * 256 + 127 * 512 + j];
  //   }
  // }
  for(int i = 0; i < g.size(); i++) { if( g[i] > 0.0001 ) ratio_gf[i] = f[i] / g[i]; }
  writeRawFile("./result/ratio_gf_float_1024-512-180.raw", ratio_gf);
  // writeRawFile("proj_n_float_512-180.raw", g);

}
