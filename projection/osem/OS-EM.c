#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 128
#define H 128
#define N 360
#define K 9	//更新回数
#define S 8 //サブセット数


const  char readFileName1[] = "backProjection_float_128-128.raw";
const  char readFileName2[] = "float_testCircleProjection_128-360_Bilinear.raw";
const  char readFileName3[] = "float_test_fproj_128-360.raw";
const  char readFileName4[] = "float_test_r_bproj_128-128.raw";


void Projecion_Bilinear_OSEM(float* f,float* g, int subset);
void backProjection_Bilinear_OSEM(float* g,float* h, int subset);
void readRawFile (const char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);
void osem(float* f, float* g, float* h);
void make_cij(float* f);


int main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * N, sizeof(float));
	float* h = (float*)calloc(H * W, sizeof(float));
	char writeFileName[50];

	readRawFile(readFileName1, sizeof(float), H * W, f);
	readRawFile(readFileName2, sizeof(float), H * N, g);

	for(int i = 0; i < K; i++)
	{
		osem(f, g, h);
		sprintf(writeFileName, "os-em%2d_float_128-128.raw", i+1);
		writeRawFile(writeFileName, sizeof(float), H * W, h);
	}
	// writeRawFile("os-em_float_128-128.raw", sizeof(float), H * W, h);

}



void Projecion_Bilinear_OSEM(float* f, float* g, int subset)
{
	for(int theta_degree = subset; theta_degree < N; theta_degree = theta_degree + S)
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

					g[W*theta_degree + j] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;		
				}
			}
		}
	}
}


void backProjection_Bilinear_OSEM(float* g, float* h, int subset)
{
	for(int theta_degree = subset; theta_degree < N; theta_degree = theta_degree + S)
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

				if(-(W - 1.0f) / 2.0f <= X && X <= (W - 1.0f) / 2.0f && -(H - 1.0f) / 2.0f <= Y && Y <= (H - 1.0f) / 2.0f )
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

void osem(float* f, float* g, float* h)
{
	float* cij_proj = (float*)calloc(H * N, sizeof(float));	
	float* cij 		= (float*)calloc(H * W, sizeof(float));	
	float* f_proj 	= (float*)calloc(H * N, sizeof(float));
	float* ratio_gf = (float*)calloc(H * N, sizeof(float));
	float* r_bproj 	= (float*)calloc(H * W, sizeof(float));

	
	//backProjection_Bilinear_OSEM(cij_proj, cij, 1);
	for(int subset = 0; subset < S; subset++)	
	{
		make_cij(cij_proj);
		backProjection_Bilinear_OSEM(cij_proj, cij, subset);
		
		Projecion_Bilinear_OSEM(f, f_proj, subset);
		/////////////////////////////////////////////////////////
		writeRawFile(readFileName3, sizeof(float), W * N, f_proj);
		/////////////////////////////////////////////////////////
		for(int j = 0; j < W * N; j++)
        {
            if(f_proj[j] > 0.0000000001)
                ratio_gf[j] =  g[j] / f_proj[j];
        }
		backProjection_Bilinear_OSEM(ratio_gf, r_bproj, subset);

		/////////////////////////////////////////////////////////////////////////////////////////

        writeRawFile(readFileName4, sizeof(float), W * H, r_bproj);

        /////////////////////////////////////////////////////////////////////////////////////////

		for(int j = 0; j < W * H; j++)
			h[j] = (r_bproj[j] * f[j])  / cij[j];
		for(int j = 0; j < W * H ; j++)
			f[j] = h[j];
	}

	free(cij_proj);
	free(cij);
	free(f_proj);
	free(ratio_gf);
	free(r_bproj);
}

void make_cij(float* f)
{
	for(int i = 0; i < W * 360; i++)
		f[i] = 1;
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
