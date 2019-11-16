#include <iostream>
#include "fileio.h"


int main()
{
  std::vector<float> f(64*64*64, 0.);
  std::vector<float> g(f.size(), 0.);
  readRawFile("Shepp_float_64-64-64.raw", f);
  for(int i = 0; i < f.size(); i++)
  {
  	if(f[i] > 15.) { g[i] = 0.15; }
  	else if(f[i] > 5.) { g[i] = 0.28; }
  	else { g[i] = 0.; }

  }
  writeRawFile("absorp_map_shepp_float_64-64-64.raw", g);

}
