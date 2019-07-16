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
  std::vector<float> f(512*256*180, 0.);

  readRawFile("primary_single_pinhole_float_512-256-180.raw", f);

	for(int i = 0; i < 180; i++)
  {
		  f[i * 512 * 256 + 0] = 0.;
	}

  writeRawFile("pprimary_single_pinhole_float_512-256-180.raw", f);
}
