#include <stdio.h>
#include <stdlib.h>
#define H 256
#define W 256


int main()
{
    unsigned char   f[H][W];       
    unsigned char   g[H][W];       
    FILE            *fp1 ,*fp2;

    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_uchar_256-256.raw" , "rb");
	fread( f , sizeof(unsigned char) , H * W , fp1);
	fclose( fp1 );
    /****************************************************************/


	for(int i = 0;i < H; i++)
	{
		for(int j = 0; j < W; j++)
		{
			if(f[i][j] < 128)
				g[i][j] = 0;
			else
				g[i][j] = 255;
		}
	}



    /**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "kadai1.raw" , "wb" );
    fwrite( g , sizeof(unsigned char) , H * W , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

