
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

Ca（0.28）
H2o（0.15）
air(0.01)

*/


void create_ellipse_h2o(float* f);

template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);



int main()
{
	float* f = (float*)calloc(H * W * D, sizeof(float));

	{
		create_ellipse_h2o(f);
		writeRawFile("mu-map_h2o_ellipse_float_128-128-128.raw", H * W * D, f);
	}

	

}

void create_ellipse_h2o(float* f)
{
	int a = 4;
	int b = 10;
	int c = 10;

	for(int i = 0; i < D; i++)
	{
		for(int j = 0; j < D; j++)
		{
			for(int k = 0; k < D; k++)
			{
				float x = - (W - 1.0) / 2.0 + k;
				float y =   (H - 1.0) / 2.0 - j;
				float z =   (D - 1.0) / 2.0 - i;
				
				// printf("x = %f, y = %f\n", x, y);
				float pixel_size = 0.25;
				x *= pixel_size;
				y *= pixel_size;
				z *= pixel_size;
				bool test = pow(x, 2.) / pow(a, 2.) + pow(y, 2.) / pow(b, 2.) + pow(z, 2.) / pow(c, 2.) < 1;
				if(test) { f[i * W * H + j * W + k] = 0.15; }
			}
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

