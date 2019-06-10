#include <stdio.h>
#include <stdlib.h>
#define N 256

void shirokuro( unsigned char f[][N], unsigned char g[][N] );

void main(void){

    int             i, j;
    unsigned char   f[N][N], g[N][N];
    FILE            *fp;


    fp = fopen( "lenna.256" , "rb" );
    if(fp == NULL){
        printf("faild to open file.");
        exit(1);
    }
    fread( f, sizeof(unsigned char), N*N, fp );
    fclose( fp );


    shirokuro(f, g);


    fp = fopen( "lenna_new.256" , "wb" );
    if(fp == NULL){
        printf("faild to open file.");
        exit(1);
    }
    fwrite( g, sizeof(unsigned char), N*N, fp );
    fclose( fp );

}

void shirokuro( unsigned char f[][N], unsigned char g[][N] ){

    int i,j;

    for(i=0; i<N; i++){
        for(j=0; j<N; j++){

            g[i][j] = 255 - f[i][j];

        }
    }

}


