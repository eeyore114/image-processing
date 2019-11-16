#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128
	
void rotateImage(float f[],float g[]);
void inputProjecion_Bilinear(float f[],float g[]/*int array[],float theta*/);


int main()
{
	float   	f[W*H];				/** 原画像用の配列 **/
	float 		g[W*360];				/** 回転後の画像の配列 **/
    int 		array[W];		    /** 投影の値を格納する配列 **/
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
	inputProjecion_Bilinear(f,g);

	/**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "CircleProjection_bilinear_float_128-360.raw" , "wb" );
    fwrite( g , sizeof(float) , 128 * 360  , fp2 );
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
				float x = - W / 2 + j;
				float y = - H / 2 - 50.0f + i;

				
				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f <= X && X <= (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f <= Y && Y <= (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 +(W - 1.0)/2.0 ;
					float I = (H - 1.0)/2.0 - y0;

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
