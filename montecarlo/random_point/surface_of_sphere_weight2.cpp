/*
	乱数生成（メルセンヌ）
	参考サイト：http://vivi.dyndns.org/tech/cpp/random.html

	0から1の乱数の書き方（uniform_real_distribution）
	http://siv3d.hateblo.jp/entry/2013/02/17/231829

*/

/*
メルセンヌ乱数を用いて半径50の球の中に点をうつ

大量に打つことにより円を生成できる


*/
#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#define N 30000000
#define H 128
#define W 128
#define L 128
#define R 50.0


template <typename T> 
void writeRawFile (const char fname[], const size_t num, T* image);

int main()
{
	float* f = (float*)calloc(H * W * L, sizeof(float));

	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<float> radius_random(0.0, 1.);        
    uniform_real_distribution<float> phi_random(0.0, 2.0 * M_PI);        
    uniform_real_distribution<float> theta_random(0.0, M_PI);        
    uniform_real_distribution<float> random(-1.0, 1.0);        
    
    for(int i = 0; i < N; i++)
    {
	    float radius = radius_random(mt);
	    float phi = phi_random(mt);
	    float theta = theta_random(mt);

	    float z = random(mt);
	    float x = sqrt(1 - z * z) * cos(phi);
	    float y = sqrt(1 - z * z) * sin(phi);

	    x *= 50;
	    y *= 50;
	    z *= 50;

	    int J = x + (W - 1.0) / 2.0;
		int I = (H - 1.0) / 2.0 - y;
		int K = (L - 1.0) / 2.0 - z;

		f[W * H * K + W * I + J] += 1.0f;

	 }

	 writeRawFile("surface_sphere_weight_float_128-128-128.raw", H * W * L, f);

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


