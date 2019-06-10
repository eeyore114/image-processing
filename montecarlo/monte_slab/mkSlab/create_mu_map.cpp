
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>

#define H 128
#define W 128
#define D 128

/*
mu-mapの作成

mu-map_float_128-128-128.raw

左半分：Ca（0.28）
右半分：H2o（0.15）

*/

void create_mu_map(float* f);
void create_mu_map_ca(float* f);
void create_mu_map_h2o(float* f);

template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);


int main()
{
	float* f = (float*)calloc(H * W * D, sizeof(float));

	create_mu_map_ca(f);

	writeRawFile("mu-map_ca_float_128-128-128.raw", H * W * D, f);

}

void create_mu_map_multi(float* f)
{
	for(int i = 0; i < D; i++)
	{
		for(int j = 0; j < D; j++)
		{
			for(int k = 0; k < D; k++)
			{
				if(k < W / 2) { f[i * W * H + j * W + k] = 0.28; }
				else		  { f[i * W * H + j * W + k] = 0.15; }
			}

		}
	}
}

void create_mu_map_ca(float* f)
{
	for(int i = 0; i < D; i++)
	{
		for(int j = 0; j < D; j++)
		{
			for(int k = 0; k < D; k++)
				f[i * W * H + j * W + k] = 0.28;
		}
	}
}

void create_mu_map_h2o(float* f)
{
	for(int i = 0; i < D; i++)
	{
		for(int j = 0; j < D; j++)
		{
			for(int k = 0; k < D; k++)
				f[i * W * H + j * W + k] = 0.15;
		}
	}
}






template <class T>
void writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}

