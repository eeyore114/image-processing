#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define W 256
#define H 256

const char readFileName[] = "lenna_uchar_256-256.raw";
const char writeFileName[] = "rotate_lenna_uchar_256-256.raw";
const float theta = 45.0f * M_PI / 180.0f;

void readRawFile (const char fname[], const size_t size, const size_t num, void* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);
void copyUchar2Float(unsigned char* src, float* dst, const size_t num);
void copyFloat2Uchar(float* src, unsigned char* dst, const size_t num);
void rotateBilinear(float* f,float* g,float theta);

int main()
{

	//原画像ファイル読み込み
	unsigned char* image = (unsigned char*)calloc(H * W, sizeof(unsigned char));
	readRawFile(readFileName,sizeof(unsigned char), H * W, image);

	//画像の配列fにimageをコピー
	//imageはいらなくなるのでfreeにしておく
	float* f = (float*)calloc(H * W, sizeof(float));
	copyUchar2Float(image, f, H * W);

	//回転後の画像を格納する配列gを準備する
	float* g = (float*)calloc(H * W, sizeof(float));

////////////////////////////////////////////////////////

	rotateBilinear(f,g,theta);

////////////////////////////////////////////////////////

	copyFloat2Uchar(g, image, H * W);

	writeRawFile(writeFileName, sizeof(unsigned char), H * W, image);

	free(image);
	free(f);
	free(g);

}

void readRawFile (const char fname[], const size_t size, const size_t num, void* image)
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

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image)
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

void copyUchar2Float(unsigned char* src, float* dst, const size_t num)
{
	for(size_t i = 0;i < num; i++){ dst[i] = (float)src[i]; }
}

void copyFloat2Uchar(float* src, unsigned char* dst, const size_t num)
{
	for(size_t i = 0;i < num; i++){ dst[i] = (unsigned char)src[i]; }
}


void rotateBilinear(float f[],float g[],float theta)
{

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			
			//回転後の画像の画素(i, j)の座標(x, y)
			float x = -(W - 1.0)/2.0 + j;
			float y = (H - 1.0)/2.0 - i;
			
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

				g[W*i+j] = f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;			
			}			
		}
	}	
}





