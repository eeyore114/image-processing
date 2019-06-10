/*
 *
 * uni_rand_wMT.c -- メルセンヌツイスタを用いて一様乱数を生成
 *
 */

#include<iostream>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define NUM 1E6
#define H 64
#define W 64
#define D 64

template <typename T> 
void writeRawFile (const char fname[], const size_t num, T* image);

int main(void){

	unsigned long seed;
	double r, theta, phi, x, y, z;
	FILE *fp;
	char filename[64];
	float* array = (float*)calloc(H * W * D, sizeof(float));


	//時定数で乱数生成器のseed値を初期化
	seed = (unsigned long)time(NULL);
	init_genrand(seed);

	//生成した乱数のヒストグラム用の出力ファイル
	snprintf(filename, sizeof(filename),  "uniform_rand_hist_3d.txt");
	fp = fopen(filename, "w");

	//単位円内に分布する乱数x, yをNUM個生成して書き出し
	for(int i = 0; i < NUM; i++){
		//条件分岐を使用しないでx, yを決定する
		// r =  genrand_real3();
		r = 1.;
	    float rand = 2 * (genrand_real3() - 0.5);
	    theta = genrand_real3() * 2 * M_PI;

		x = cbrt(r) * sqrt(1 - rand * rand) * cos(theta);
		y = cbrt(r) * sqrt(1 - rand * rand) * sin(theta);
		z = cbrt(r) * rand;

		x /= 0.05;
		y /= 0.05;
		z /= 0.05;
		fprintf(fp, "%lf %lf %lf\n", x, y, z);

		int J = x + (W - 1.0) / 2.0;
		int I = (H - 1.0) / 2.0 - y;
		int K = (D - 1.0) / 2.0 - z;

		array[K * W * H + I * W + J]++;
	}
	writeRawFile("rand_sphere_float_64-64-64.raw", W * H * D, array);

	fclose(fp);

	printf("%s is generated!\n", filename);

	return 0;
}

template <typename T> 
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

