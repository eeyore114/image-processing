
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>
#define SIZE_MUMAP 128



template <class T>
void writeRawFile (string fname, const size_t num, T* image);

int main()
{
  float* f = (float*)calloc(SIZE_MUMAP * SIZE_MUMAP * SIZE_MUMAP, sizeof(float));

  for(int i = 0; i < SIZE_MUMAP * SIZE_MUMAP * SIZE_MUMAP; i++) { f[i] = 0.; }

  writeRawFile("mu_map_air_float_128-128-128.raw", SIZE_MUMAP * SIZE_MUMAP * SIZE_MUMAP, f);
}


template <class T>
void writeRawFile(string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname.c_str());
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}
