/*
	乱数生成（メルセンヌ）
	参考サイト：http://vivi.dyndns.org/tech/cpp/random.html

	0から1の乱数の書き方（uniform_real_distribution）
	http://siv3d.hateblo.jp/entry/2013/02/17/231829

*/


#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#define N 10000000
#define SIZE 128


template <typename T> 
void writeRawFile (const char fname[], const size_t num, T* image);


int main()
{
	float pt = 0;
	float* f = (float*)calloc(SIZE * SIZE, sizeof(float));


	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<float> radius_random(0.0,2500.0); //0からr^2までの範囲       
    uniform_real_distribution<float> theta_random(0.0, 2.0 * M_PI);        
    

    //ランダムに点を打つ
    for(int i = 0; i < N; i++)
    {

    	float radius = radius_random(mt);
	    float theta = theta_random(mt);

	    float x = sqrt(radius) * cos(theta);
	    float y = sqrt(radius) * sin(theta);


	    int J = x + (SIZE - 1.0) / 2.0;
		int I = (SIZE - 1.0) / 2.0 - y;

		f[I * SIZE + J] += 1.;
	}

	writeRawFile("circle_random_wight_float_128-128.raw", SIZE * SIZE, f);
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




