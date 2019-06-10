
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <random>
void mkSurfaceSourceImg(float* f, int surface_source_w, int surface_source_h);


template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);


int main()
{
	int surface_source_w = 2048;
	int surface_source_h = 1024;
	float* surface_source = (float*)calloc(surface_source_w * surface_source_w, sizeof(float));

	mkSurfaceSourceImg(surface_source, surface_source_w, surface_source_h);

	writeRawFile("surfaceSource_float_2048-1024.raw", surface_source_w * surface_source_w, surface_source);
}

void mkSurfaceSourceImg(float* f, int surface_source_w, int surface_source_h)
{
	for(int i = 0; i < surface_source_h; i++)
	{
		for(int j = 0; j < surface_source_w; j++)
		{
			float x = - (surface_source_w - 1.) / 2. + j;
			float y = (surface_source_h - 1.) / 2. - i;
			int r = 10;
			r /= 0.5;
			if(abs(x) > r || abs(y) > r) { continue; }
			f[i * surface_source_w + j] = 1;
		}
	}
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
