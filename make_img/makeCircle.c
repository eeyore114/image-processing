#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void makeCircleImage(float* f);

int main()
{
	int W = 128;
	int H = 128;
	float   	f[W * H] = {};
	FILE    	*fp;
	makeCircleImage(f);

	fp = fopen( "circle_float_128-128.raw" , "wb" );
    fwrite( f , sizeof(float) , 128 * 128  , fp );
    fclose( fp );

    return 0;
}

void makeCircleImage(float* f)
{
	int W = 128;
	int H = 128;
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{	
			//(i, j)座標を(x, y)座標に変換
			float x = -(W - 1.0) / 2.0 + j;
			float y = (H - 1.0) / 2.0 - i;

			//円の中なら100の画素値を入れる
			if(pow(x, 2.0f) + pow(y, 2.0f) < pow(50, 2))
				f[W * i + j] = 100.0f;
			else
				f[W * i + j] = 0.0f;
		}
	}	
}

