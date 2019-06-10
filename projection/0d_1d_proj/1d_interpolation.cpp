#include <iostream>
using namespace std;
#include <stdio.h>
#include <math.h>


void interpolation_1d(float* f, float* g);


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image);

template <typename T> void 
writeRawFile (const char fname[], const size_t num, T* image);


int main()
{
	int H = 128;
	int W = 128;
	int N = 360;

	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	readRawFile("circle.raw", H * W, f);

	interpolation_1d(f, g);

	writeRawFile("1d_proj_float_128-360.raw", H * N, g);




}

void interpolation_1d(float* f, float* g)
{
	int H = 128;
	int W = 128;
	int N = 360;

	for(int theta_degree = 0; theta_degree < N; theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for(int i = 0; i < H; i++)
		{
			for(int j = 0; j < W; j++)
			{
				float x = - (W - 1.0f) / 2.0f + j;
				float y =   (H - 1.0f) / 2.0f - i;

				float s =   x * cosf(theta) + y * sinf(theta);
	            float t = - x * sinf(theta) + y * cosf(theta);

	            //ij軸に変換して小数点を切り捨てする
	            float J = s + (W - 1) / 2;
	            int idx = floor(J);

				//j + 1番目の検出器までの距離をdxとする。
				float dx = J - idx; 

    			g[W * theta_degree + idx]	  += fabsf((1-dx)) * f[W * i + j];
    			g[W * theta_degree + idx + 1] += dx  * f[W * i + j];
        		
			}
		}
	}
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




template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image)
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
