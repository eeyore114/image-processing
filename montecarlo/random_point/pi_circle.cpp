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


int main()
{
	float pt = 0;

	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<float> score(0.0,1.0);        // [0.0, 1.0] 範囲の一様乱数
    
    for(int i = 0; i < N; i++)
    {
	    float x = score(mt);
	    float y = score(mt);
	    if((pow(x, 2.0f) + pow(y, 2.0f)) <= 1.)
	    	pt++;
	 }

	 float pi = 4.0f * pt / N;


	 cout << "case : number of points =  " << N << endl;
	 cout << "pi = " << pi << endl;
}

