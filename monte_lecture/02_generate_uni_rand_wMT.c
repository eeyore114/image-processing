/*
 *
 * generate_uni_rand_wMT.c -- メルセンヌツイスタを用いて一様乱数を生成
 *
 */

#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define NUM 1E5

int main(void){

	unsigned long seed;
	double r, x, y;
	FILE *fp;
	char filename[64];

	//時定数で乱数生成器のseed値を初期化
	seed = (unsigned long)time(NULL);
	init_genrand(seed);

	//生成した乱数のヒストグラム用のcsvファイル
	snprintf(filename, sizeof(filename),  "uniform_rand_hist.txt");
	fp = fopen(filename, "w");
	int n = 0;

	//区間(-1,1)の乱数x, yをNUM個生成して書き出し
	for(int i = 0; i < NUM; i++){
//		r = genrand_real3();
		x = 2 * ( genrand_real3() - 0.5 ); //区間(-1, 1)の乱数xを生成
		y =	2 * ( genrand_real3() - 0.5 ); //区間(-1, 1)の乱数yを生成
		if( x * x + y * y < 1 ) {				//単位円内に入っている場合
			n++;
			fprintf(fp, "%lf %lf\n", x, y);
		}
	}
	printf("pi = %f\n", 4 * n / NUM);

	fclose(fp);

	printf("%s is generated!\n", filename);

	return 0;
}

