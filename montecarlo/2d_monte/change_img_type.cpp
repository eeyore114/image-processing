#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>


template <class T>
void readRawFile (std::string fname, const size_t num, T* image);
template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);


int main()
{
	double* f = (double*)calloc(180 * 180, sizeof(double));
	float* g = (float*)calloc(180 * 180, sizeof(float));

	readRawFile("cylinder_h2o_single_pinhole_00_double_180-180.raw", 180 * 180, f);

	for(int i = 0; i < 180 * 180; i++) { g[i] = (float)f[i]; }

	writeRawFile("cylinder_h2o_single_pinhole_00_float_180-180.raw", 180 * 180, g);

}


template <class T>
void readRawFile (std::string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname.c_str());
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}

template <class T>
void writeRawFile(std::string fname, const size_t num, T* image)
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
