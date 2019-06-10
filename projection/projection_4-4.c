#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 4
#define H 4
#define theta_degree 90.0f

const float theta = theta_degree * M_PI / 180.0f;
void inputProjecion(float f[],int array[],float theta);

int main()
{
	float   	f[W*H]={1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};       /** 原画像用の配列 **/
    int 		array[W];		  /** 投影の値を格納する配列 **/
    FILE    	*fp1 ,*fp2;
    	
    /**原画像ファイル読み込み********************************/
	/*fp1 = fopen("gazo_4-4.raw" , "rb");
	fread( f , sizeof(float) , W * H , fp1);
	fclose( fp1 );*/
	
	//arrayの初期化
	for (int i = 0; i < H; ++i)
	{
		array[i] = 0;
	}

	//fの値確認
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			if(j == W-1)
				printf("	f[%d]	= %.0f\n", W*i + j,f[W*i + j]);
			else
				printf("	f[%d]	= %.0f", W*i + j,f[W*i + j]);
		}
	}

	printf("%.0f\n\n",theta_degree);

	inputProjecion(f,array,theta);

	//pの値表示
	for (int i = 0; i < H; i++)
	{
		printf("array[%d] = %d\n",i, array[i]);	
	}

	/**トリミング画像ファイル書き込み**************************************/
   /* fp2 = fopen( "float_kadai2_45.raw" , "wb" );
    fwrite( g , sizeof(float) , 256 * 256  , fp2 );
    fclose( fp2 );*/
    /****************************************************************/

    return 0;
}

void inputProjecion(float f[],int array[],float theta)
{

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{	
			//(i, j)座標を(x, y)座標に変換
			float x = -(W - 1.0)/2.0 + j;
			float y = (H - 1.0)/2.0 - i;

			//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
			float S = x * cosf(theta) - y * sinf(theta);
			float T = x * sinf(theta) + y * cosf(theta);

			int index = floor(S + W / 2);
			array[index] += f[W*i + j];
		}
	}	
}



