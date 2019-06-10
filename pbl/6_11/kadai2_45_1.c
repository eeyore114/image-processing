#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 256
#define H 256
#define zoom_i 0.2//i軸方向の倍率
#define zoom_j 0.2//j軸方向の倍率

const float theta = 45.0f * M_PI / 180.0f;
void rotateBilinear(float f[],float g[]);

int main()
{
	float   f[W*H];       /** 原画像用の配列 **/
    float   g[W*H];       /** トリミング後の画像の配列 **/
    FILE    *fp1 ,*fp2;
    	
    /**原画像ファイル読み込み********************************/
	fp1 = fopen("lenna_float_256-256.raw" , "rb");
	fread( f , sizeof(float) , 256 * 256 , fp1);
	fclose( fp1 );

	for (int i = 0; i < H; ++i)
	{
		for (int j = 0; j < W; ++j)
		{
			g[W*i+j] = 0;
		}
	}

	rotateBilinear(f,g);

	/**トリミング画像ファイル書き込み**************************************/
    fp2 = fopen( "float_kadai2_45.raw" , "wb" );
    fwrite( g , sizeof(float) , 256 * 256  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;



}

void rotateBilinear(float f[],float g[])
{

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			
			//回転後の画像の画素(i, j)の座標(x, y)
			float x = -(W - 1.0)/2.0 + j;
			float y = (H - 1.0)/2.0 - i;
			
			//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
			float X = x / zoom_i;
			float Y = y / zoom_j;

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

				g[W*i+j] = f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;			
			}			
		}
	}	
}





