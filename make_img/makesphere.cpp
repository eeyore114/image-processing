#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void makeCircleImage(float* f);

int main()
{
	int W = 128;
	int H = 128;
	int L = 128;
	float* f = (float*)calloc(H * W * L, sizeof(float));
	FILE    	*fp;
	makeCircleImage(f);

	fp = fopen( "sphere_float_128-128-128.raw" , "wb" );
    fwrite( f , sizeof(float) , H * W * L  , fp );
    fclose( fp );

    return 0;
}

void makeCircleImage(float* f)
{
	int W = 128;
	int H = 128;
	int L = 128;
	float voxel_size = 0.2;
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			for(int k = 0; k < L; k++)
			{	
				//(i, j)座標を(x, y)座標に変換
				float x = - (W - 1.0) / 2.0 + j;
				float y =   (H - 1.0) / 2.0 - i;
				float z = - (L - 1.0) / 2.0 + k;

				x *= voxel_size;
				y *= voxel_size;
				z *= voxel_size;

				//円の中なら100の画素値を入れる
				if(pow(x, 2.0f) + pow(y, 2.0f) + pow(z, 2.0f) < pow(10, 2) )
					f[H * W * k + W * i + j] = 100.0f;
				else
					f[H * W * k +  W * i + j] = 0.0f;
			}
		}
	}	
}

