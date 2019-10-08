#include <iostream>
#include <vector>
#include "fileio.h"
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)

int main()
{
  int w = 512;
  int h = 256;
  int d = 180;
  int after_d = 3;
  std::vector<float> f(w * h * d, 0.);
  std::vector<float> g(w * h * after_d, 0.);
  readRawFile("./img/primary_detector_float_512-256-180.raw", f);
  int delta = d / after_d;
  rep(k, after_d)rep(i, h)rep(j, w) g[k * w * h + i * w + j] = f[(delta * k) * w * h + i * w + j];
  writeRawFile("./result/primary_detector_float_512-256-3.raw", g);
}
