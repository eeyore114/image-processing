#include <iostream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 512
#define H 512
#define N 360
#define k 100	//更新回数

char readFileName[] = "Brain_float_512-512.raw";


void Projecion_Bilinear(float f[],float g[]);
void backProjection_Bilinear(float g[],float h[]);
void mlem(float* f, float* g, float* h);
void make_cij(float* f);


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image);

template <typename T> void 
writeRawFile (const char fname[], const size_t num, T* image);


int main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * N, sizeof(float));
	float* h = (float*)calloc(H * W, sizeof(float));
	char writeFileName[50];

	readRawFile(readFileName, H * W, f);

	Projecion_Bilinear(f, g);


	for(int i = 0; i < H * W; i++)
		f[i] = 0.;

	backProjection_Bilinear(g, f);

	for(int i = 0; i < k; i++)
	{
		cout << i << endl;
		mlem(f, g, h);
		// sprintf(writeFileName, "ML-EM%02d.raw", i+1);
		// writeRawFile(writeFileName, sizeof(float), H * W, h);
		if(i == 30 - 1 ) { writeRawFile("ml-em_30_float_512-512.raw", H * W, h); }
		if(i == 40 - 1 ) { writeRawFile("ml-em_40_float_512-512.raw", H * W, h); }
		if(i == 60 - 1 ) { writeRawFile("ml-em_60_float_512-512.raw", H * W, h); }
		if(i == 70 - 1 ) { writeRawFile("ml-em_70_float_512-512.raw", H * W, h); }
		if(i == 80 - 1 ) { writeRawFile("ml-em_80_float_512-512.raw", H * W, h); }
		if(i == 90 - 1 ) { writeRawFile("ml-em_90_float_512-512.raw", H * W, h); }
		// if(i == 100 - 1 ) { writeRawFile("ml-em_100_float_512-512.raw", H * W, h); }
	}

}



void mlem(float* f, float* g, float* h)
{
	float* cij_proj = (float*)calloc(H * 360, sizeof(float));	
	float* cij 		= (float*)calloc(H * W, sizeof(float));	
	float* f_proj 	= (float*)calloc(H * 360, sizeof(float));
	float* ratio_gf = (float*)calloc(H * 360, sizeof(float));
	float* r_bproj 	= (float*)calloc(H * W, sizeof(float));

	make_cij(cij_proj);
	backProjection_Bilinear(cij_proj, cij);
	Projecion_Bilinear(f, f_proj);

	for(int i = 0; i < W * 360; i++)
		ratio_gf[i] =  g[i] / f_proj[i] ;

	backProjection_Bilinear(ratio_gf, r_bproj);

	for(int i = 0; i < W * H; i++)
		h[i] = (r_bproj[i] * f[i])  / cij[i];

	for(int i = 0; i < W * H ; i++)
		f[i] = h[i];
}

void make_cij(float* f)
{
	for(int i = 0; i < W * 360; i++)
		f[i] = 1.;
}



void Projecion_Bilinear(float f[],float g[])
{
	for(int theta_degree = 0;theta_degree < 360;theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50; i++)
		{
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - W / 2 + j;
				float y = - H / 2 - 50.0f + i;


				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 +(W - 1.0)/2.0 ;
					float I = (H - 1.0)/2.0 - y0;

					//indexは配列の番号
					int index1 = W * I + J;
					int index2 = W * I + J + 1;
					int index3 = W * ( I + 1) + J;
					int index4 = W * ( I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);

					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下

					g[W*theta_degree + j] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;		
				}
			}
		}
	}
}


void backProjection_Bilinear(float g[],float h[])
{


	for(int theta_degree = 0;theta_degree < 360;theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50; i++)
		{
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - W / 2 + j;
				float y = - H / 2 - 30.0f + i;

				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 + (W - 1.0) / 2.0 ;
					float I = (H - 1.0) / 2.0 - y0;

					//indexは配列の番号
					int index1 = W * I + J;
					int index2 = W * I + J + 1;
					int index3 = W * (I + 1) + J;
					int index4 = W * (I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);

					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下

					h[index1] += g[W * theta_degree + j] * S4;
					h[index2] += g[W * theta_degree + j] * S3;
					h[index3] += g[W * theta_degree + j] * S2;
					h[index4] += g[W * theta_degree + j] * S1;
				}
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
