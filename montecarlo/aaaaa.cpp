#include <iostream>
#include <stdlib.h>
#include <math.h>

template <class T>
void readRawFile (std::string fname, const size_t num, T* image);

template <class T>
void writeRawFile(std::string fname, const size_t num, T* image);

int main()
{
	float* f = (float*)calloc(180*180*180, sizeof(float));
	float* g = (float*)calloc(180*180*180, sizeof(float));
	readRawFile("primary_float_180-180-180.raw", 180*180*180, f);

	for(int i = 0; i < 180; i++)
	{
		for(int j = 0; j < 180; j++)
		{
			for(int k = 0; k < 180; k++)
			{
				g[k * 180 * 180 + i * 180 + j] = f[0 * 180 * 180 + i * 180 + j];
			}
		}
	}

	writeRawFile("primary2_float_180-180-180.raw", 180*180*180, g);

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
