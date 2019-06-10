#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128
#define L 128
#define N 360
#define D 512

const char readFileName[] = "sphere_float_128-128-128.raw";
const char writeFileName1[] = "sphere_singlepinhole_float_512-256-360.raw";
const char writeFileName2[] = "sphere_singlepinhole_float_128-128-128.raw";

void Projecion_SinglePinhole_3d(float* f,float* g);
void backProjection_Singlepinhole_3d(float* f, float* g);
void search_LessThan30_3d(float* f, float* g);


template <typename T>
void readRawFile (const char fname[], const size_t num, T* image);

template <typename T>
void writeRawFile (const char fname[], const size_t num, T* image);

int main(void)
{
	float* f = (float*)calloc(H * W * L, sizeof(float));
	
	float* g = (float*)calloc(512 * 256 * 360, sizeof(float));
	float* h = (float*)calloc(H * W * L, sizeof(float));

	readRawFile(readFileName, H * W * L, f);

	Projecion_SinglePinhole_3d(f, g);
	writeRawFile(writeFileName1,  512 * 256 * 360, g);

	backProjection_Singlepinhole_3d(g, h);
	writeRawFile(writeFileName2,  128 * 128 * 128, h);
}


void Projecion_SinglePinhole_3d(float* f,float* g)
{
	//ピクセルサイズは0.17 cm × 0.17 cm
	float pixel_size = 0.17f;
	float* theta_collimator_xy = (float*)calloc( 512 * 256, sizeof(float));
	float* theta_collimator_zy = (float*)calloc( 512 * 256, sizeof(float));

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy);

	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
	{
		cout << theta_degree << "° : projection processing..." << endl;
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = - 25.0f;

		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		float z = 0.;

		for(int m = 0; m < 256; m++)
		{
			for(int n = 0; n < 512; n++)
			{
				//スケーリング処理
				float X = x / pixel_size;
				float Y = y / pixel_size;
				float Z = z / pixel_size;

				if(theta_collimator_xy[m * 512 + n] != 0.0f && theta_collimator_zy[m * 512 + n] != 0.0f)
				{
					for (int sp = 0; sp < 300; sp++)
					{
						if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f &&
						   -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f &&
						   -(L - 1.0f) / 2.0f < Z && Z < (L - 1.0f) / 2.0f )
						{
							//双線形補完を行う左上の画素の座標(x0,y0)
							float x0 = floor(X - 0.5f) + 0.5f;
							float y0 = ceil(Y - 0.5f) + 0.5f;
							float z0 = ceil(Z - 0.5f) + 0.5f;

							//i,jに戻す
							float J = x0 + (W - 1.0f) / 2.0f;
							float I = (H - 1.0f) / 2.0f - y0;
							float K = (L - 1.0f) / 2.0f - z0;

							int index1 = W * H * K + W * I + J;
							int index2 = W * H * K + W * I + (J + 1);
							int index3 = W * H * K + W * (I + 1) + J;
							int index4 = W * H * K + W * (I + 1) + (J + 1);
							int index5 = W * H * (K + 1) + W * I + J;
							int index6 = W * H * (K + 1) + W * I + (J + 1);
							int index7 = W * H * (K + 1) + W * (I + 1) + J;
							int index8 = W * H * (K + 1) + W * (I + 1) + (J + 1);

							//(X , Y)から左上までの距離dx, dy, dz
							float dx = fabs(X - x0);
							float dy = fabs(Y - y0);
							float dz = fabs(Z - z0);

							float V1 = dx * dy * dz;
							float V2 = (1.0f - dx) * dy * dz;
							float V3 = dx * (1.0f - dy) * dz;
							float V4 = (1.0f - dx) * (1.0f - dy) * dz;
							float V5 = dx * dy * (1.0f - dz);
							float V6 = (1.0f - dx) * dy * (1.0f - dz);
							float V7 = dx * (1.0f - dy) * (1.0f - dz);
							float V8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);

							g[512 * 256 * theta_degree + 512 * m + n] += f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 +
							 									         f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;

						}

						X +=  cos(theta_collimator_zy[512 * m + n]) * sin(theta - theta_collimator_xy[512 * m + n]);
						Y +=  cos(theta_collimator_zy[512 * m + n]) * cos(theta - theta_collimator_xy[512 * m + n]);
						Z +=  -sin(theta_collimator_zy[512 * m + n]);
					}
				}
			}
		}
	}
}



void backProjection_Singlepinhole_3d(float* f, float* g)
{
	//ピクセルサイズは0.17 cm × 0.17 cm
	float pixel_size = 0.17f;
	float* theta_collimator_xy = (float*)calloc( 512 * 256, sizeof(float));
	float* theta_collimator_zy = (float*)calloc( 512 * 256, sizeof(float));

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy);

	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
	{
		cout << theta_degree << "° : backprojection processing..." << endl;
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -25.0f;

		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		float z = 0.;

		for(int m = 0; m < 256; m++)
		{
			for(int n = 0; n < 512; n++)
			{
				//スケーリング処理
				float X = x / pixel_size;
				float Y = y / pixel_size;
				float Z = z / pixel_size;


				if(theta_collimator_xy[m * 512 + n] != 0.0f && theta_collimator_zy[m * 512 + n] != 0.0f)
				{
					for (int sp = 0; sp < 300; sp++)
					{
						if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f &&
						   -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f &&
						   -(L - 1.0f) / 2.0f < Z && Z < (L - 1.0f) / 2.0f )
						{
							//双線形補完を行う左上の画素の座標(x0,y0)
							float x0 = floor(X - 0.5f) + 0.5f;
							float y0 = ceil(Y - 0.5f) + 0.5f;
							float z0 = ceil(Z - 0.5f) + 0.5f;

							//i,jに戻す
							float J = x0 + (W - 1.0f) / 2.0f;
							float I = (H - 1.0f) / 2.0f - y0;
							float K = (L - 1.0f) / 2.0f - z0;

							int index1 = W * H * K + W * I + J;
							int index2 = W * H * K + W * I + (J + 1);
							int index3 = W * H * K + W * (I + 1) + J;
							int index4 = W * H * K + W * (I + 1) + (J + 1);
							int index5 = W * H * (K + 1) + W * I + J;
							int index6 = W * H * (K + 1) + W * I + (J + 1);
							int index7 = W * H * (K + 1) + W * (I + 1) + J;
							int index8 = W * H * (K + 1) + W * (I + 1) + (J + 1);

							//(X , Y)から左上までの距離dx, dy, dz
							float dx = fabs(X - x0);
							float dy = fabs(Y - y0);
							float dz = fabs(Z - z0);

							float V1 = dx * dy * dz;
							float V2 = (1.0f - dx) * dy * dz;
							float V3 = dx * (1.0f - dy) * dz;
							float V4 = (1.0f - dx) * (1.0f - dy) * dz;
							float V5 = dx * dy * (1.0f - dz);
							float V6 = (1.0f - dx) * dy * (1.0f - dz);
							float V7 = dx * (1.0f - dy) * (1.0f - dz);
							float V8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);


							g[index1] += f[512 * 256 * theta_degree + 512 * m + n] * V8;
							g[index2] += f[512 * 256 * theta_degree + 512 * m + n] * V7;
							g[index3] += f[512 * 256 * theta_degree + 512 * m + n] * V6;
							g[index4] += f[512 * 256 * theta_degree + 512 * m + n] * V5;
							g[index5] += f[512 * 256 * theta_degree + 512 * m + n] * V4;
							g[index6] += f[512 * 256 * theta_degree + 512 * m + n] * V3;
							g[index7] += f[512 * 256 * theta_degree + 512 * m + n] * V2;
							g[index8] += f[512 * 256 * theta_degree + 512 * m + n] * V1;
						}

						X +=  cos(theta_collimator_zy[512 * m + n]) * sin(theta - theta_collimator_xy[512 * m + n]);
						Y +=  cos(theta_collimator_zy[512 * m + n]) * cos(theta - theta_collimator_xy[512 * m + n]);
						Z +=  -sin(theta_collimator_zy[512 * m + n]);
					}
				}
			}
		}
	}
}

void search_LessThan30_3d(float* f, float* g)
{
	//確認用
	float* h = (float*)calloc(512 * 256, sizeof(float));

	for(int i = 0; i < 256; i++)
	{
		for(int j = 0; j < 512; j++)
		{
			float i0 = i * 0.08f + 0.04f;
			float j0 = j * 0.08f + 0.04f;
			float x = - 512 * 0.08f / 2 + j0;
			float z =   256 * 0.08f / 2 - i0;
			float theta_Ditector_xy = atanf(x / 7.5f);
			float theta_Ditector_zy = atanf(z / 7.5f);

			// 半径7.5/sqrt(3)はピンホールに対する角度が30の時の円の半径
			float radius = 7.5 / sqrt(3);

			if(pow(x, 2.0f) + pow(z, 2.0f) < pow(radius, 2.0f))
			{
				f[i * 512 + j] = theta_Ditector_xy;
				g[i * 512 + j] = theta_Ditector_zy;
				h[i * 512 + j] = 100.0f;//確認用
			}
			else
			{
				f[i * 512 + j] = 0.0f;
				g[i * 512 + j] = 0.0f;
				h[i * 512 + j] = 0.;//確認用
			}
		}
	}
	writeRawFile("circleLessThan30_float_512-256.raw", 512 * 256, h);
}


template <typename T>
void readRawFile (const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}


template <typename T>
void writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}
