#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define H 	128
#define W 	128
#define PI 3.14

unsigned char readfilename[] 	= "circle.raw";
unsigned char powerspectrum[] 	= "powerspectrum.raw";
unsigned char writefilename[] 	= "DFT_circle.raw";

void  DFT(unsigned char fname1[], unsigned char fname2[]);
void IDFT(unsigned char fname1[], unsigned char fname2[]);
void readRawFile (const unsigned char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const unsigned char fname[], const size_t size, const size_t num, float* image);

void main(void)
{
	DFT(readfilename, powerspectrum);
//	IDFT(powerspectrum, writefilename);
}	

void DFT(unsigned char fname1[], unsigned char fname2[])
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* R = (float*)calloc(H * W, sizeof(float));
	float* I = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	readRawFile(fname1, sizeof(float), H * W, f);
	//FIXME : xyuvの座標の表現
	for(int i = 0; i < W; i++)
	{
		for (int k = 0; k < W; k++)
		{
			int u = k - W / 2;
			for(int j = 0; j < H; j++)
		 	{
				for(int l = 0; l < H; l++)
				{
					int v = l - W / 2;
					R[i * H + j] +=  f[k * H + l] * cos(2 * PI * (u * j / W + v * i / H));
					I[i * H + j] += -f[k * H + l] * sin(2 * PI * (u * j / W + v * i / H));
				}	
			}
		}
	}

	for (int i = 0; i < H * W; i++)
		g[i] = pow(R[i], 2.0f) + pow(I[i], 2.0f);
	
	writeRawFile(fname2, sizeof(float), H * W, g);

	free(f);
	free(R);
	free(I);
	free(g);
}

void IDFT(unsigned char fname1[], unsigned char fname2[])
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* R = (float*)calloc(H * W, sizeof(float));
	float* I = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	readRawFile(fname1, sizeof(float), H * W, f);
	for(int i = 0; i < W; i++)
	{
		for (int k = 0; k < W; k++)
		{
			int u = k - W / 2;
			for(int j = 0; j < H; j++)
		 	{
				for(int l = 0; l < H; l++)
				{
					int v = l - W / 2;
					R[k * H + l] += f[i * H + j] * cos(2 * PI * (u * j / W + v * i / H));
					I[k * H + l] += f[i * H + j] * sin(2 * PI * (u * j / W + v * i / H));
				}	
			}
		}
	}

	for (int i = 0; i < H * W; i++)
	{
		g[i] = R[i] / (H * W);
	}
	writeRawFile(fname2, sizeof(float), H * W, g);
}

void readRawFile (const unsigned char fname[], const size_t size, const size_t num, float* image)
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

void writeRawFile(const unsigned char fname[], const size_t size, const size_t num, float* image)
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
