/*
 *
 * uni_rand_wMT.c -- メルセンヌツイスタを用いて一様乱数を生成
 *
 */

#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define NUM 1E3

int main(void){

	unsigned long seed;
	double r, theta, x, y;
	FILE *fp;
	char filename[64];

	//時定数で乱数生成器のseed値を初期化
	seed = (unsigned long)time(NULL);
	init_genrand(seed);

	//生成した乱数のヒストグラム用の出力ファイル
	snprintf(filename, sizeof(filename),  "uniform_rand_hist.txt");
	fp = fopen(filename, "w");

	//単位円内に分布する乱数x, yをNUM個生成して書き出し
	for(int i = 0; i < NUM; i++){
		//条件分岐を使用しないでx, yを決定する
		r = genrand_real3();
	    theta = genrand_real3() * 2 * M_PI;
		// x = sqrt(r) * cos(theta);
		x = r * cos(theta);
		// y = sqrt(r) * sin(theta);
		y = r * sin(theta);
		fprintf(fp, "%lf %lf\n", x, y);
	}

	fclose(fp);

	printf("%s is generated!\n", filename);

	return 0;
}

