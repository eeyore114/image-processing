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
#define SIZE 128

long long N = 99999999;

template <typename T> 
void writeRawFile (const char fname[], const size_t num, T* image);


int main()
{
	float pt = 0;
	float* f = (float*)calloc(SIZE * SIZE * SIZE, sizeof(float));


	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<float> score(-64.0,64.0);        // [0.0, 1.0] 範囲の一様乱数
    

    //ランダムに点を打つ
    for(int i = 0; i < N; i++)
    {
	    float x = score(mt);
	    float y = score(mt);
	    float z = score(mt);

	    int J = x + (SIZE - 1.0) / 2.0;
		int I = (SIZE - 1.0) / 2.0 - y;
		int K = (SIZE - 1.0) / 2.0 - z;

		f[SIZE * SIZE * K + SIZE * I + J] += 1.;
	}


    for(int i = 0; i < SIZE; i++)
    {
    	for(int j = 0; j < SIZE; j++)
    	{
    		for(int k = 0; k < SIZE; k++)
    		{
				float x = j - (SIZE - 1.0) / 2.0;
				float y = (SIZE - 1.0) / 2.0 - i;
				float z = (SIZE - 1.0) / 2.0 - k;

				if(pow(x, 2.) + pow(y, 2.) + pow(z, 2.) > pow(50, 2.))
				{
					int J = x + (SIZE - 1.0) / 2.0;
					int I = (SIZE - 1.0) / 2.0 - y;
					int K = (SIZE - 1.0) / 2.0 - z;
					f[SIZE * SIZE * K + SIZE * I + J] = 0.;
				}
				cout << SIZE * SIZE * k + SIZE * i + j << " : processing..." << endl;;
    		}
    	}
    }


	writeRawFile("sphere_random_float_128-128-128.raw", SIZE * SIZE * SIZE, f);
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




