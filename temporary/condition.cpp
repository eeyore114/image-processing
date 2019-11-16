#include <iostream>
#include "fileio.h"

int main()
{
	std::vector<float> v(128*128*128, 0);
	readRawFile("sphere_ML-EM100_float_128-128-128.raw", v);
	for(int i = 0; i < v.size(); i++) { if(v[i] > 50) { v[i] = 0; } }
	writeRawFile("resphere_ML-EM100_float_128-128-128.raw", v);
}