#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128

char readFileName[] = "circle_float_128-128.raw";
char writeFileName[] = "backProjection_float_128-128.raw";

void readRawFile (char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(char fname[], const size_t size, const size_t num, float* image);
void Projection_Bilinear(float f[],float g[]);
void backProjection_Bilinear(float f[],float g[]);


int main()
{

	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * 360, sizeof(float));
	float* h = (float*)calloc(H * W, sizeof(float));
	readRawFile(readFileName,sizeof(float), H * W, f);



////////////////////////////////////////////////////////
	
	Projection_Bilinear(f,g);
	writeRawFile("projection_float_128-360.raw", sizeof(float), H * 360, g);
	backProjection_Bilinear(g,h);

////////////////////////////////////////////////////////


	writeRawFile(writeFileName, sizeof(float), H * W, h);

	free(f);
	free(g);
	free(h);


    return 0;
}

void Projection_Bilinear(float f[],float g[])
{
	

	for(int theta_degree = 0;theta_degree < 360; theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50 ; i++)
		{ 
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - W / 2 + j;
				float y = - H / 2 - 30.0f + i;

				
				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 + ( W - 1.0 )/2.0 ;
					float I = ( H - 1.0 ) / 2.0 - y0;

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

					g[W * theta_degree + j] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;		
				}
			}
		}
	}
}

void backProjection_Bilinear(float g[],float h[])
{
	

	for(int theta_degree = 0;theta_degree < 360; theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < H + 50; i++)
		{ 
			for (int j = 0; j < W; j++)
			{
				//検出器の座標を定義
				float x = - W / 2 + j;
				float y = - H / 2 - 30.0f + i;

				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 + (W - 1.0) / 2.0 ;
					float I = (H - 1.0) / 2.0 - y0;

					//indexは配列の番号
					int index1 = W * I + J;
					int index2 = W * I + J + 1;
					int index3 = W * (I + 1) + J;
					int index4 = W * (I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);
				
					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下

					h[index1] += g[W * theta_degree + j] * S4;
					h[index2] += g[W * theta_degree + j] * S3;
					h[index3] += g[W * theta_degree + j] * S2;
					h[index4] += g[W * theta_degree + j] * S1;
				}
			}
		}
	}
}


void readRawFile (char fname[], const size_t size, const size_t num, float* image)
{
	//ファイルを開く
	FILE* fp = fopen(fname,"rb");

	//ファイルを開くことができなかった場合のエラー処理
	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	//データの読み込み
	size_t ret = fread(image, size, num,fp);

	//データを読み込むことができなかった場合のエラー処理
	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);

	}

	//ファイルを閉じる
	fclose(fp);
}

void writeRawFile(char fname[], const size_t size, const size_t num, float* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, size, num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);


}

