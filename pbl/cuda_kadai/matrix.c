#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

const int N = 1024;

void matrix_cal(float* a, float* b, float* c);

int main()
{
    // 領域の確保
    float* a = (float*)calloc(N * N, sizeof(float));
    float* b = (float*)calloc(N * N, sizeof(float));
    float* c = (float*)calloc(N * N, sizeof(float));
    
    // 数値の初期化
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            
            a[N * i + j] = i;
            b[N * i + j] = j;
            c[N * i + j] = 0.;
            
        }
    }

    matrix_cal(a, b, c);

 /*   for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            printf("c[%d] = %f\n", i * N + j, c[i * N + j]);
    }
    */
    

}


void matrix_cal(float* a, float* b, float* c)
{
    struct timeval startTime, endTime;  // 構造体宣言
    clock_t startClock, endClock;       // clock_t型変数宣言

    gettimeofday(&startTime, NULL);     // 開始時刻取得
    startClock = clock();               // 開始時刻のcpu時間取得



///////////////////////////////////////////////////////////////////

    // 行列の計算
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++){
                c[N * i + j] += a[N * i + k] * b[N * k + j];
            }
        }
    }

//////////////////////////////////////////////////////////////////

    gettimeofday(&endTime, NULL);       // 開始時刻取得
    endClock = clock();                 // 開始時刻のcpu時間取得
    
    time_t diffsec = difftime(endTime.tv_sec, startTime.tv_sec);    // 秒数の差分を計算
    suseconds_t diffsub = endTime.tv_usec - startTime.tv_usec;      // マイクロ秒部分の差分を計算

     double realsec = (diffsec+diffsub*1e-6) * 1000;                          // 実時間を計算

     printf("time = %lf[ms]\n", realsec);

}
