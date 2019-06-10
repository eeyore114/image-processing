#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define H 	128
#define W 	128
#define PI  3.14159265358979


char readfilename[] 	= "circle.raw";
char powerspectrum[] 	= "powerspectrum3.raw";
char writefilename[] 	= "DFT_circle3.raw";

void DFT(char fname1[], char fname2[], float* Re, float* Im);
void IDFT(float* Re, float* Im, char fname2[]);
void readRawFile ( char fname[],  size_t size,  size_t num, float* image);
void writeRawFile( char fname[],  size_t size,  size_t num, float* image);


int main(void)
{
	float* Re = (float*)calloc(H * W, sizeof(float));
	float* Im = (float*)calloc(H * W, sizeof(float));

	DFT(readfilename, powerspectrum, Re, Im);
	IDFT(Re, Im, writefilename);

	free(Re);
	free(Im);
}

void DFT(char fname1[], char fname2[], float* Re, float* Im)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	readRawFile(fname1, sizeof(float), H * W, f);

	for(int i = 0; i < H; i++)
	{
		int v = - i + H / 2;
		for(int j = 0; j < W; j++)
		{
			int u = j - W / 2;
			for(int k = 0; k < H; k++)
			{
				int y = - k + H / 2;
				for(int m = 0; m < W; m++)
				{
					int x = m - W / 2;
					Re[i * W + j] +=  f[k * W + m] * cos(2 * PI * (u * x / W + v * y / H));
					Im[i * W + j] += -f[k * W + m] * sin(2 * PI * (u * x / W + v * y / H));
				}
			}
			g[i * W + j] += log10(pow(Re[i * W + j], 2.0f) + pow(Im[i * W + j], 2.0f));
		}
	}

	writeRawFile(fname2, sizeof(float), H * W, g);

	free(f);
	free(g);
}


void IDFT(float* Re, float* Im, char fname2[])
{
	float* g = (float*)calloc(H * W, sizeof(float));

	for(int i = 0; i < H; i++)
	{
		int v =  - i + H / 2;
		for(int j = 0; j < W; j++)
		{
			int u = j - W / 2;
			for(int k = 0; k < H; k++)
			{
				int y = - k + H / 2;
				for(int m = 0; m < W; m++)
				{
					int x = m - W / 2;
					g[i * W + j] +=  Re[k * W + m] * cos(2 * PI *(u * x / W + v * y / H)) / (W * H) - Im[k * W + m] * sin(2 * PI * (u * x / W + v * y / H)) / (W * H);
				}
			}
		}
	}
	writeRawFile(fname2, sizeof(float), H * W, g);

	free(g);
}



void readRawFile (char fname[], size_t size,  size_t num, float* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, size, num, fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}


void writeRawFile( char fname[],  size_t size,  size_t num, float* image)
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