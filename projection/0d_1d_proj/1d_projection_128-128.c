#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128


void Projecion_NearestNeighbor(float f[],float g[];/*int array[],float theta*/);

int main()
{
	float   	f[W*H];				/** 原画像用の配列 **/
	float 		g[W*360];			/** 投影を行った値を格納する配列 **/
    int 		array[W];		  /** 投影の値を格納する配列 **/
    FILE    	*fp1 ,*fp2;
    	
    /**原画像ファイル読み込み********************************/
	fp1 = fopen("circle_float_128-128.raw" , "rb");
	fread( f , sizeof(float) , W * H , fp1);
	fclose( fp1 );
	
	//arrayの初期化
	for (int i = 0; i < H; ++i)
	{
		array[i] = 0;
	}

	for(int i = 0; i < W * 360; i++)
	{
		g[i] = 0;
	}

	Projecion_NearestNeighbor(f,g/*,theta_degree*/);


	/**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "float_testCircleProjection_128-360_1d.raw" , "wb" );
    fwrite( g , sizeof(float) , 128 * 360  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

void Projecion_NearestNeighbor(float f[],float g[])
{

	for(int theta_degree = 0;theta_degree < 360;theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < W; j++)
			{	
				//(i, j)座標を(x, y)座標に変換
				float x = -(W - 1.0)/2.0 + j;
				float y = (H - 1.0)/2.0 - i;

				//座標(x, y)を原点を中心に-theta回転させた座標(S, T)
				float S = x * cosf(theta) - y * sinf(theta);
				float T = x * sinf(theta) + y * cosf(theta);

				
				int index = floor(S + W / 2);
				float ds = index + 0.5 - (S + W / 2);
				if(ds > 0)
				{
					g[W*theta_degree + index - 1] += f[W*i + j] * ds;
					g[W*theta_degree + index] 	  += f[W*i + j] * (1.0f - ds);
				}
				else
				{
					ds *= -1.0f;
					g[W*theta_degree + index] 	  += f[W*i + j] * (1.0f - ds);
					g[W*theta_degree + index + 1] += f[W*i + j] * ds;
				}

			}
		}
	}	
}



