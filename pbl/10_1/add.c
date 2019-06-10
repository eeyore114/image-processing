#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#define N 1024*255

  void main()
{
	struct timeval startTime, endTime;  // 構造体宣言
    clock_t startClock, endClock;       // clock_t型変数宣言

    gettimeofday(&startTime, NULL);     // 開始時刻取得
    startClock = clock();               // 開始時刻のcpu時間取得

    /////////////////////////////////////////////////////////////////////

	int a[N] = {};

	for(int i = 0; i < N; i++)
		a[i] += i; 

	/////////////////////////////////////////////////////////////////////

	gettimeofday(&endTime, NULL);       // 開始時刻取得
    endClock = clock();                 // 開始時刻のcpu時間取得
    
    time_t diffsec = difftime(endTime.tv_sec, startTime.tv_sec);    // 秒数の差分を計算
    suseconds_t diffsub = endTime.tv_usec - startTime.tv_usec;      // マイクロ秒部分の差分を計算

	 double realsec = diffsec+diffsub*1e-6;                          // 実時間を計算

	 printf("time = %lf[s]\n", realsec);

}
