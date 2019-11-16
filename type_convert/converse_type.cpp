#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#define W 128
#define H 128

const  char readFileName[] = "shepp_double_128-128.raw";
const  char writeFileName[] = "shepp_float_128-128.raw";


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image);

void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);

int main()
{
	float* g = (float*)calloc(H * W, sizeof(float));
	{
		double* f = (double*)calloc(H * W, sizeof(double));
	
		readRawFile(readFileName, H * W, f);	
		for(int i = 0; i < W * H; i++)
		{
			g[i] = (float)f[i];
		}
	}

	writeRawFile(writeFileName, sizeof(float), H * W, g);

	
}


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}


void writeRawFile(const char fname[], const size_t size, const size_t num, float* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, size, num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}
