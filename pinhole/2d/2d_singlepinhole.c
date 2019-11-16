/*長さはあってる*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128
#define N 360
#define D 512

const char readFileName[] = "circle.raw";
const char writeFileName1[] = "float_singlepinhole_512-360.raw";
const char writeFileName2[] = "float_singlepinhole_128-128.raw";

void Projecion_Bilinear_SinglePinhole(float* f,float* g);
void backProjection_Bilinear_Singlepinhole(float* f, float* g);
void search_LessThan30(float* f);
void readRawFile (const char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);




int main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(D * N, sizeof(float));
	float* h = (float*)calloc(H * W, sizeof(float));

	readRawFile(readFileName,sizeof(float), H * W, f);

	Projecion_Bilinear_SinglePinhole(f, g);
	writeRawFile(writeFileName1, sizeof(float), D * N, g);

	backProjection_Bilinear_Singlepinhole(g, h);
	writeRawFile(writeFileName2, sizeof(float), H * W, h);
}


void Projecion_Bilinear_SinglePinhole(float* f,float* g)
{
	//ピクセルサイズは0.17 cm × 0.17 cm
	float pixel_size = 0.17f;
	float* h = (float*)calloc( D, sizeof(float));

	search_LessThan30(h);

	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		
		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -25.0f;

		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta); 

		for(int m = 0; m < D; m++)
		{
			//スケーリング処理
			float X = x / pixel_size;
			float Y = y / pixel_size;
			if(h[m] != 0.0f)
			{
				for (int n = 0; n < H + W; n++)
				{
					if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
					{
						//双線形補完を行う左上の画素の座標(x0,y0)
						float x0 = floor(X - 0.5) + 0.5;
						float y0 = ceil(Y - 0.5) + 0.5;

						//i,jに戻す
						float J = x0 + (W - 1.0)/2.0;
						float I = (H - 1.0)/2.0 - y0;

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

						g[D * theta_degree + m] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;		
					}
					
					/*
					これを変形したのが下の式
					X +=  cos(M_PI / 2 - (theta - h[m]);
					Y +=  sin(M_PI / 2 - (theta - h[m]); */

					X +=  sin(theta - h[m]);
					Y +=  cos(theta - h[m]); 

					if(m == 256 && n == H + W - 1 && theta_degree == 0)
						printf("X : %f 	Y : %f\n", X, Y);	
				}
			}
			else
			{
				g[D * theta_degree + m] = 0.0f;
			}	
				
		}
	}
}



void backProjection_Bilinear_Singlepinhole(float* f, float* g)
{
	//ピクセルサイズは0.17 cm × 0.17 cm
	float pixel_size = 0.17f;
	float* h = (float*)calloc( D, sizeof(float));

	search_LessThan30(h);

	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		
		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -25.0f;

		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta); 

		for(int m = 0; m < D; m++)
		{
			//スケーリング処理
			float X = x / pixel_size;
			float Y = y / pixel_size;

			if(h[m] != 0.0f)
			{
				for (int n = 0; n < H + W; n++)
				{
					if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f )
					{
						//双線形補完を行う左上の画素の座標(x0,y0)
						float x0 = floor(X - 0.5) + 0.5;
						float y0 = ceil(Y - 0.5) + 0.5;

						//i,jに戻す
						float J = x0 + (W - 1.0)/2.0;
						float I = (H - 1.0)/2.0 - y0;

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

						g[index1] += f[D * theta_degree + m] * S4;
						g[index2] += f[D * theta_degree + m] * S3;
						g[index3] += f[D * theta_degree + m] * S2;
						g[index4] += f[D * theta_degree + m] * S1;	
					}
					
					/*
					これを変形したのが下の式
					X +=  cos(M_PI / 2 - (theta - h[m]);
					Y +=  sin(M_PI / 2 - (theta - h[m]);
					 */

					X +=  sin(theta - h[m]);
					Y +=  cos(theta - h[m]); 
				}
			}
		}
	}
}



void search_LessThan30(float* f)
{
		int count = 0;

	for(int j = 0; j < D; j++)
	{
		float j0 = j * 0.08f + 0.04f;
		float x = - D * 0.08f / 2 + j0;
		float theta_Ditector = atanf(x / 7.5f);
		
		if(fabsf(theta_Ditector) < M_PI / 6)
		{
			count++;
			f[j] = theta_Ditector;
		}
		else
		{
			f[j] = 0.0f;
		}
	}
}




void readRawFile (const char fname[], const size_t size, const size_t num, float* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, size, num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}


void writeRawFile(const char fname[], const size_t size, const size_t num, float* image)
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


