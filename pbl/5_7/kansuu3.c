#include <stdio.h>
#include <stdlib.h>

void shirokuro(unsigned char *f,unsigned char *g){
	int i,j;
	
	for(i=0;i<256;i++){
		for(j=0;j<256;j++){
			g[i*256 + j] = 255 - f[i*256 + j];
		}
	}
	
}



int main()
{
    int             i, j;
    unsigned char   *f;       /** 原画像用の配列 **/
    unsigned char   *g;       /** トリミング後の画像の配列 **/
    FILE            *fp1 ,*fp2;

    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_uchar_256-256.raw" , "rb");
	fread( f , 1 , 256 * 256 , fp1);
	fclose( fp1 );

    /****************************************************************/

	shirokuro(f,g);


    /**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "kannsuu3.raw" , "wb" );
    fwrite( g , 1 , 256 * 256  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

