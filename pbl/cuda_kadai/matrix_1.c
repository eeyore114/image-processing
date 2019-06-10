#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 256;

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

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            printf("c[%d] = %f\n", i * N + j, c[i * N + j]);
    }
    
    

}


void matrix_cal(float* a, float* b, float* c)
{
    // 行列の計算
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++){
                c[N * i + j] += a[N * i + k] * b[N * k + j];
            }
        }
    }

}
