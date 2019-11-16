#include <iostream>
#include <vector>
#include "fileio.h"
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)

int main()
{
  int w = 512;
  int h = 256;
  int d = 180;
  int after_d = 90;
  std::vector<float> f1(w * h * d, 0.);
  std::vector<float> f2(w * h * d, 0.);
  std::vector<float> g1(w * h * after_d, 0.);
  std::vector<float> g2(w * h * after_d, 0.);
  readRawFile("./img/primary_detector_float_512-256-180.raw", f1);
  readRawFile("./img/primary_detector_rectangle_float_512-256-180.raw", f2);
  int delta = d / after_d;
  rep(k, after_d)rep(i, h)rep(j, w) g1[k * w * h + i * w + j] = f1[(delta * k) * w * h + i * w + j];
  rep(k, after_d)rep(i, h)rep(j, w) g2[k * w * h + i * w + j] = f2[(delta * k) * w * h + i * w + j];
  writeRawFile("./result/primary_detector_float_512-256-90.raw", g1);
  writeRawFile("./result/primary_detector_rectangle_float_512-256-90.raw", g2);
}
