#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>

#define SIZE 128

// void scale_by2(float** f, float** g);

template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image);

template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image);


int main()
{

   

    float f[SIZE][SIZE] = {};
    float g[2 * SIZE][2 * SIZE] = {};
    FILE            *fp1, *fp2;



    // readRawFile("Brain_float_128-128.raw", SIZE * SIZE, f);

    fp1 = fopen( "Brain_float_128-128.raw" , "rb" );
    fread( f , sizeof(float) , SIZE * SIZE , fp1 );
    fclose( fp1 );

    // scale_by2(f, g);

    for(int i = 0; i < SIZE; i++)
    {
        for(int j = 0; j < SIZE; j++)
        {
            g[i * 2][j * 2] = f[i][j];
            g[i * 2][j * 2 + 1] = (f[i][j] + f[i][j + 1]) / 2;
            g[i * 2 + 1][j * 2] = (f[i][j] + f[i + 1][j]) / 2;
            g[i * 2 + 1][j * 2 + 1] = (f[i][j] + f[i][j + 1] + f[i + 1][j] + f[i + 1][j + 1]) / 4;
        }
    }

    // writeRawFile("Brain_float_256-256.raw", 2 * SIZE * 2 * SIZE, g);

    fp2 = fopen( "Brain_float_256-256.raw" , "wb" );
    fwrite( g , sizeof(float) , 4 * SIZE * SIZE , fp2 );
    fclose( fp2 );

    return 0;
}

/*
void scale_by2(float** f, float** g)
{
    for(int i = 0; i < SIZE; i++)
    {
        for(int j = 0; j < SIZE; j++)
        {
            g[i * 2][j * 2] = f[i][j];
            g[i * 2][j * 2 + 1] = (f[i][j] + f[i][j + 1]) / 2;
            g[i * 2 + 1][j * 2] = (f[i][j] + f[i + 1][j]) / 2;
            g[i * 2 + 1][j * 2 + 1] = (f[i][j] + f[i][j + 1] + f[i + 1][j] + f[i + 1][j + 1]) / 4;
        }
    }
}
*/




template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image)
{
    FILE* fp = fopen(fname,"rb");

    if(fp == NULL)
    {
        printf("failed to open %s.\n",fname);
        exit(-1);
    }

    size_t ret = fread(image, sizeof(T), num,fp);

    if(num != ret)
    {
        printf("failed to read %s\n", fname);
        fclose(fp);
        exit(-1);
    }

    fclose(fp);
}


template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image)
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
        printf("failed to write %s.\n", fname);
        fclose(fp);
        exit(-1);
    }

    fclose(fp);
}

