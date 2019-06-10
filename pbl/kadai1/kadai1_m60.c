#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 256
#define H 256

void rotateCoordinate(unsigned char f[],unsigned char g[],float theta){

	int i, j;
	int	I, J;
    float x;
	float y;
	float X;
	float Y;
	float rad = theta * M_PI / 180.0;



	/***i,jをx,yに変換***/
	for (i = 0; i < H; i++)
	{
		for (j = 0; j < W; j++)
		{
			x = -(W - 1.0)/2.0 + j;
			y = (H - 1.0)/2.0 - i;
			X = x * cosf(rad) - y * sinf(rad);
			Y = x * sinf(rad) + y * cosf(rad);
			I = (H - 1.0)/2.0 - Y;
			J = (W - 1.0)/2.0 + X;
			if(I <= 256 && I >= 0 && J >= 0 && J < 256)
			g[W*I+J] = f[W*i+j];

		}
	}
	
	
}




int main()
{

	int             i, j;
	unsigned char   f[256*256];       /** 原画像用の配列 **/
    unsigned char   g[256*256];       /** トリミング後の画像の配列 **/
    FILE            *fp1 ,*fp2;
    float theta;
	theta = -60.0;
	
    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_uchar_256-256.raw" , "rb");
	fread( f , 1 , 256 * 256 , fp1);
	fclose( fp1 );

	for (i = 0; i < H; ++i)
	{
		for (j = 0; j < W; ++j)
		{
			g[W*i+j] = 0;
		}
	}

	rotateCoordinate(f,g,theta);

	/**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "kadai1_m60.raw" , "wb" );
    fwrite( g , 1 , 256 * 256  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;



}


