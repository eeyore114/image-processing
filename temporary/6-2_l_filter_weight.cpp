/******************************************************************************
        6-2_l_filter_weight.c
******************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <vector>
#include "fileio.h"

int main()
{
  int w = 512;
  int h = 256;
  std::vector<float> f(w * h, 0.);
  std::vector<float> g(f.size(), 0.);

  readRawFile("efficiency_map_float_512-256.raw", f);

	for(int i = 1; i < h - 1; i++)
  {
		for(int j = 1; j < w - 1; j++)
    {
			g[i * w + j] = (f[(i-1) * w + (j-1)] + f[(i-1) + w + j] * 2
                      + f[(i-1) * w + (j+1)] + f[i * w + j-1] * 2
                      + f[i * w + j] * 4
                      + f[i * w + j + 1] * 2
                      + f[(i + 1) * w + (j - 1)]
                      + f[(i + 1) * w + j] * 2
                      + f[(i+1) * w + j + 1]) / 16;

    }
	}

  writeRawFile("map_float_512-256.raw", g);
}

void smoothing_filter_with_weight(std::vector<float> &f, std::vector<float> &g, Condition cond)
{
	for(int i = 1; i < cond.detector_size_h - 1; i++)
  {
		for(int j = 1; j < cond.detector_size_w - 1; j++)
    {
			g[i * w + j] = (f[(i-1) * w + (j-1)] + f[(i-1) + w + j] * 2
                      + f[(i-1) * w + (j+1)] + f[i * w + j-1] * 2
                      + f[i * w + j] * 4
                      + f[i * w + j + 1] * 2
                      + f[(i + 1) * w + (j - 1)]
                      + f[(i + 1) * w + j] * 2
                      + f[(i+1) * w + j + 1]) / 16;

    }
	}
}
