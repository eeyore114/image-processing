#include <iostream>
#include "fileio.h"


int main()
{
  int w = 512;
  int h = 256;
  int d = 180;
  std::vector<float> f(w * h * d, 0.);
  std::vector<int> g(w * h, 0.);
  std::vector<int> result(g.size(), -1);
  readRawFile("./img/primary_single_pinhole_float_512-256-180.raw", f);
  readRawFile("./img/fov_11pinhole_square_5mm_int_512-256.raw", g);
  for(int i = 0; i < g.size(); i++) {	if(g[i] == -1 && f[41 * w * h + i] > 0.1) { result[i] = 1; } }
  writeRawFile("./result/fail_int_512-256.raw", result);

}
