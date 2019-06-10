#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 256
#define H 256

void LocalBinaryPattern(float f[],float g[]);

int main()
{
	float   f[W*H];       					/** 原画像用の配列 **/
    float   g[W*H];	
    FILE    *fp1 ,*fp2;
    	
    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_float_256-256.raw" , "rb");
	fread( f , sizeof(float) , 256 * 256 , fp1);
	fclose( fp1 );

	for (int i = 0; i < H ; ++i)
	{
		for (int j = 0; j < W ; ++j)
		{
			g[W*i+j] = 0;
		}
	}

	LocalBinaryPattern(f,g);

	/**トリミング画像ファイル書き込み*************************************/
    fp2 = fopen( "float_LBP.raw" , "wb" );
    fwrite( g , sizeof(float) , W*H  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

void LocalBinaryPattern(float f[],float g[])
{
	int n = 3;//3×3の画素領域を抽出
	
	for (int i = 1; i < H - 1; i++)
	{
		for (int j = 1; j < W - 1; j++)
		{
			float h[n * n];

			for(int l = 0;l < n;l++)
			{
				h[l] = f[W * (i-1) + (j - 1 + l)];
			}


			h[0] = f[W * (i-1) + (j-1)];
			h[1] = f[W * (i-1) + (j)];
			h[2] = f[W * (i-1) + (j+1)];
			h[3] = f[W * (i) + (j+1)];
			h[4] = f[W * (i+1) + (j+1)];
			h[5] = f[W * (i+1) + (j)];
			h[6] = f[W * (i+1) + (j-1)];
			h[7] = f[W * (i) + (j-1)];
			
			for(int k = n - 2; k >= 0 ; k--)
			{
				if(h[(n - 2)-k] >= f[W * i + j])
					g[W*i + j] += pow(2,k);				
			}				
		}
	}
}





