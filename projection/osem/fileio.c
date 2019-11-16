#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128

void readRawFile (const char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);

int main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	const char readFileName[] = "circle_float_128-128.raw";
	const char writeFileName[] = "result_float_128-128.raw";

	readRawFile(readFileName, sizeof(float), H * W, f);

	for(int i = 0; i < W * H; i++) { g[i] = f[i]; }

	writeRawFile(writeFileName, sizeof(float), H * W, g);
}

void readRawFile (const char fname[], const size_t size, const size_t num, float* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, size, num,fp);

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
