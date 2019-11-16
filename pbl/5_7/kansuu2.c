#include <stdio.h>
#include <stdlib.h>
#define H 256
#define W 256

void shirokuro(unsigned char f[H][W],unsigned char g[H][W]);



int main()
{
    unsigned char   f[H][W];       
    unsigned char   g[H][W];       
    FILE            *fp1 ,*fp2;

    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_uchar_256-256.raw" , "rb");
	if(fp1 == NULL)
	{
		printf("failed to open file.");
		exit(1);
	}
	fread( f , sizeof(unsigned char) , H * W , fp1);
	fclose( fp1 );
    /****************************************************************/

	shirokuro(f, g);

    /**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "kannsuu2.raw" , "wb" );
    if(fp2 == NULL)
    {
    	printf("failed to open file.");
    	exit(1);
    }
    fwrite( g , sizeof(unsigned char) , H * W  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}


void shirokuro(unsigned char f[H][W],unsigned char g[H][W])
{
	for(int i = 0; i < H; i++)
	{
		for(int j = 0;j < W; j++)
			g[i][j] = 255 - f[i][j];
	}
	
}
