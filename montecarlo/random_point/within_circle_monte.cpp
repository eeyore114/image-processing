/*
	乱数生成（メルセンヌ）
	参考サイト：http://vivi.dyndns.org/tech/cpp/random.html

	0から1の乱数の書き方（uniform_real_distribution）
	http://siv3d.hateblo.jp/entry/2013/02/17/231829

*/

/*
メルセンヌ乱数を用いて半径50の円の中に点をうつ

大量に打つことにより円を生成できる


*/
#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#define N 300000
#define H 128
#define W 128

template <typename T> 
void writeRawFile (const char fname[], const size_t num, T* image);

int main()
{
	float* f = (float*)calloc(H * W, sizeof(float));

	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<float> radius_random(0.0, 50.0);        
    uniform_real_distribution<float> theta_random(0.0, 2.0 * M_PI);        
    
    for(int i = 0; i < N; i++)
    {
	    float radius = radius_random(mt);
	    float theta = theta_random(mt);

	    float x = radius * cos(theta);
	    float y = radius * sin(theta);

	    int J = x + (W - 1.0) / 2.0;
		int I = (H - 1.0) / 2.0 - y;

		f[W * I + J] += 1.0f;

	 }


	 writeRawFile("circle_float_128-128.raw", H * W, f);

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


