
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

*/


void create_mu_map_multi(float* f);
void create_mu_map_h2o(float* f);

template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);



int main()
{
	float* f = (float*)calloc(H * W * D, sizeof(float));
	float* g = (float*)calloc(H * W * D, sizeof(float));

	// create_mu_map_h2o(f);
	create_mu_map_multi(g);

	// writeRawFile("mu-map_h2o_sphere_float_128-128-128.raw", H * W * D, f);
	writeRawFile("mu-map_multi_sphere_float_128-128-128.raw", H * W * D, g);
}

void create_mu_map_multi(float* f)
{
	float ratio = 0.2;

	for(int i = 0; i < D; i++)
	{
		for(int j = 0; j < H; j++)
		{
			for(int k = 0; k < W; k++)
			{
				float x = - (W - 1.0) / 2.0 + k;
				float y =   (H - 1.0) / 2.0 - j;
				float z =   (D - 1.0) / 2.0 - i;
				
				// printf("x = %f, y = %f\n", x, y);
				x *= ratio;
				y *= ratio;
				z *= ratio;
				if(pow(x, 2.) + pow(y, 2.) + pow(z, 2.) < pow(10, 2.)) { f[i * W * H + j * W + k] = 0.15; }
			}
		}
	}

	{

		float x = 2.;
		float y = 3.;
		float z = 3.;

		int r = 3;

		int small_size = 2 * r / ratio;

		int index1 = (H - 1.0) / 2.0 - y / ratio;
		int index2 = (W - 1.0) / 2.0 + x / ratio;
		int index3 = (H - 1.0) / 2.0 - z / ratio;

		// cout << "index1 = " << index1;
		// cout << " index2 = " << index2 << endl;

		for(int i = index3; i < index3 + small_size; i++)
		{
			for(int j = index1; j < index1 + small_size; j++)
			{
				for(int k = index2; k < index2 + small_size; k++)
				{

					
					float small_x = - (small_size - 1.0) / 2.0 + (k - index2);
					float small_y =   (small_size - 1.0) / 2.0 - (j - index1);
					float small_z =   (small_size - 1.0) / 2.0 - (i - index3);

					cout << "small_x = " << small_x;
					cout << " small_y = " << small_y << endl;

					if(pow(small_x, 2.) + pow(small_y, 2.) + pow(small_z, 2.) < pow(r / ratio, 2.)) { f[i * W * H + j * W + k] = 0.28; }
				}
			}
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
			{
				float x = - (W - 1.0) / 2.0 + k;
				float y =   (H - 1.0) / 2.0 - j;
				float z =   (D - 1.0) / 2.0 - i;
				
				// printf("x = %f, y = %f\n", x, y);
				x *= 0.2;
				y *= 0.2;
				z *= 0.2;

				if(pow(x, 2.) + pow(y, 2.) + pow(z, 2.) < pow(10, 2.)) { f[i * W * H + j * W + k] = 0.15; }
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

