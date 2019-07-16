#include <stdio.h>
#include <stdlib.h>

#define SIZE 128

int main()
{

    int             i, j;
		float temp;
    // float   f[ SIZE *SIZE * SIZE ];/*入力画像*/
    //  float   g[ 64* (SIZE / 2) * (SIZE / 2) ];/*縮小画像*/
		float* f = (float*)calloc(128*128*128, sizeof(float));
		float* g = (float*)calloc(64*64*64, sizeof(float));

    FILE            *fp1, *fp2;

    /**トリミング画像ファイル読み込み********************************/
    fp1 = fopen( "Shepp_float_128-128-128.raw" , "rb" );
    fread( f , sizeof(float) , SIZE * SIZE * SIZE , fp1 );
    fclose( fp1 );
    /****************************************************************/

		float conf = 0.;
		for(int i = 0; i < 128*128*128; i++)
		{
			if(f[i] > conf) { conf = f[i]; }
		}
		if(conf < 0.01) { printf("all 0\n"); }


    /**↓↓プログラム作成開始↓↓************************************/
		for(int d = 0; d < 64 ; d++)
		{
			for(i = 0;i<64;i++)
			{
        for(j = 0;j<64;j++)
				{
					temp = f[d * 2 * 128 * 128 + i*2 * 128 + j*2];
					if(temp < f[d * 2 * 128 * 128 + i*2 * 128 + j*2+1])
					temp = f[d * 2 * 128 * 128 + i*2 * 128 + j*2+1];
					if(temp < f[d * 2 * 128 * 128 + (i*2+1) * 128 + j*2])
					temp = f[d * 2 * 128 * 128 + (i*2+1) * 128 + j*2];
					if(temp < f[d * 2 * 128 * 128 + (i*2+1) * 128 + j*2+1])
					temp = f[d * 2 * 128 * 128 + (i*2+1) * 128 + j*2+1];

          g[d * 64*64+ i* 64 + j] = temp;
					if(d * 64*64+ i* 64 + j > 131072 - 100 && d * 64*64+ i* 64 + j < 131072+ 10)
					 printf("%f\n", temp);

        }
    	}
		}




    /**↑↑プログラム作成終了↑↑************************************/


    /**縮小画像ファイル書き込み**************************************/
	fp2 = fopen( "Shepp_float_64-64-64.raw" , "wb" );
    fwrite( g , sizeof(float) ,  64*64*64, fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;

}
