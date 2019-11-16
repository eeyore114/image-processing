#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
using namespace std;

int height = 512;
int width = 512;

void Projection_0d_interpolation(float* f, float* g);
void backProjection_0d_interpolation(float* f, float* g);

template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image);

template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image);

int main()
{
    int H = height;
    int W = width;
    int N = 360;
	float* f = (float*)calloc(H * W, sizeof(float));
    float* g = (float*)calloc(H * N, sizeof(float));
    float* h = (float*)calloc(H * W, sizeof(float));
	readRawFile("Brain_float_512-512.raw", H * W, f);	

	Projection_0d_interpolation(f, g);

	writeRawFile("0d_proj_float_512-360.raw", H * N, g);

    backProjection_0d_interpolation(g, h);

    writeRawFile("0d_bproj_float_512-512.raw", H * W, h);
}


void Projection_0d_interpolation(float* f, float* g)
{
    int H = height;    //heiht
    int W = width;    //width
    int N = 360;    //theta(degree) 0 ~ 360

    for(int theta_degree = 0; theta_degree < 360; theta_degree++)
    {
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H; i++)
		{
    		for (int j = 0; j < W; j++)
        	{
        		//画像の画素(i, j)の座標(x, y)
        		float x =  - ( W - 1.0f ) / 2.0f + j;
       			float y =    ( H - 1.0f ) / 2.0f - i;
        
        		// 座標(x, y)を原点中心に-theta回転させた座標(X, Y)
        		float s =   x * cosf(theta) + y * sinf(theta);
        		float t = - x * sinf(theta) + y * cosf(theta);

        		int idx = floor(s + W / 2);
        		
        		g[W * theta_degree + idx] += f[W * i + j];
        	}
		}
	}
}


void backProjection_0d_interpolation(float* f, float* g)
{
    int H = height;    //heiht
    int W = width;     //width
    int N = 360;       //theta(degree) 0 ~ 360

    for(int theta_degree = 0; theta_degree < 360; theta_degree++)
    {
        const float theta = theta_degree * M_PI / 180.0f;
        for (int i = 0; i < H; i++)
        {
            for (int j = 0; j < W; j++)
            {
                //画像の画素(i, j)の座標(x, y)
                float x =  - ( W - 1.0f ) / 2.0f + j;
                float y =    ( H - 1.0f ) / 2.0f - i;
        
                // 座標(x, y)を原点中心に-theta回転させた座標(X, Y)
                float s =   x * cosf(theta) + y * sinf(theta);
                float t = - x * sinf(theta) + y * cosf(theta);

                int idx = floor(s + W / 2);
                
                g[W * i + j] += f[W * theta_degree + idx];
            }
        }
    }
}


template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s.\n",fname);
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
		printf("failed to write %s.\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}
