#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128

void inputProjecion_Bilinear(float f[],float g[]);
void backProjection_Bilinear(float g[],float h[]);

int main()
{
	float   	f[W*H] ={};				/** 原画像用の配列 **/
	float 		g[W*360] = {};			/** 投影を行った値を格納する配列 **/
	float   	h[W*H] = {};				/** 原画像用の配列 **/
    FILE    	*fp1 ,*fp2;
    	
    /**原画像ファイル読み込み********************************/
	fp1 = fopen("circle_float_128-128.raw" , "rb");
	fread( f , sizeof(float) , W * H , fp1);
	fclose( fp1 );

	inputProjecion_Bilinear(f,g);
	backProjection_Bilinear(g,h);

    /**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "backProjection_float_128-128.raw" , "wb" );
    fwrite( h , sizeof(float) , W*H  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

void inputProjecion_Bilinear(float f[],float g[])
{
	for(int theta_degree = 0;theta_degree < 360;theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50; i++)
		{ 
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - (W / 2) + j;
				float y = - (H / 2) - 30.0f + i;

				
				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f <= X && X <= (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f <= Y && Y <= (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floorf(X - 0.5) + 0.5;
					float y0 = ceilf(Y - 0.5) + 0.5;

					//i,jに戻す
					int J = x0 +(W - 1.0)/2.0 ;
					int I = (H - 1.0)/2.0 - y0;

					//indexは配列の番号
					int index1 = W * I + J;
					int index2 = W * I + J + 1;
					int index3 = W * ( I + 1) + J;
					int index4 = W * ( I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);
				
					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下

					g[W*theta_degree + j] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;		
				}
			}
		}
	}
}

void backProjection_Bilinear(float g[],float h[])
{
	for(int theta_degree = 0;theta_degree < 360;theta_degree++)
	{
		float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50; i++)
		{ 
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - (W / 2.0f) + j;
				float y = - (H / 2.0f) - 30.0f + i;

				
				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f < X && X <= (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f <= Y && Y < (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floorf(X - 0.5f) + 0.5f;
					float y0 = ceilf(Y - 0.5f) + 0.5f;

					//i,jに戻す
					int J = (W - 1.0)/2.0 + x0;
					int I = (H - 1.0)/2.0 - y0;

					//indexは配列の番号
					int index1 = W * I + J;
					int index2 = W * I + J + 1;
					int index3 = W * ( I + 1 ) + J;
					int index4 = W * ( I + 1 ) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);
				
					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下


					
					h[index1] += S4 * g[W*theta_degree + j];
					h[index2] += S3 * g[W*theta_degree + j];
					h[index3] += S2 * g[W*theta_degree + j];
					h[index4] += S1 * g[W*theta_degree + j];

				}
			}
		}
	}
}


