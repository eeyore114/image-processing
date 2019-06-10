#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define H 	128
#define W 	128
#define PI 3.1415

unsigned char readfilename[] 	= "circle.raw";
unsigned char powerspectrum[] 	= "powerspectrum.raw";
unsigned char writefilename[] 	= "DFT_circle.raw";

void DFT(unsigned char fname1[], unsigned char fname2[], float* Re, float* Im);
void dft_1d_w(float* f, float* g, float* Re, float* Im, int i);
void dft_1d_h(float* f, float* g, float* Re, float* Im, int j);
void IDFT(float* Re, float* Im, unsigned char fname2[]);
void idft_1d_w(float* Re, float* Im, float* g, int i);
void idft_1d_h(float* Re, float* Im, float* g, int j);
void idft_1d_w2(float* Re, float* Im, float* g, int i);
void idft_1d_h2(float* Re, float* Im, float* g, int j);

void readRawFile (const unsigned char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const unsigned char fname[], const size_t size, const size_t num, float* image);


void main(void)
{
	float* Re = (float*)calloc(H * W, sizeof(float));
	float* Im = (float*)calloc(H * W, sizeof(float));

	DFT(readfilename, powerspectrum, Re, Im);
	IDFT(Re, Im, writefilename);

	free(Re);
	free(Im);
}

void DFT(unsigned char fname1[], unsigned char fname2[], float* Re, float* Im)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	readRawFile(fname1, sizeof(float), H * W, f);

	for(int i = 0; i < H; i++)
	{
		dft_1d_w(f, g, Re, Im, i);
	}

	for(int j = 0; j < W; j++)
	{
		dft_1d_h(f, g, Re, Im, j);
	}

	writeRawFile(fname2, sizeof(float), H * W, g);

	free(f);
	free(g);
}

void dft_1d_w(float* f, float* g, float* Re, float* Im, int i)
{
	float* re = (float*)calloc(W, sizeof(float));
	float* im = (float*)calloc(W, sizeof(float));

	for(int j = 0; j < W; j++)
	{
		int u = j - W / 2;
		for(int x = 0; x < W; x++)
		{
			re[j] +=  f[i * W + x] * cos(2 * M_PI * u * x / W);
			im[j] += -f[i * W + x] * sin(2 * M_PI * u * x / W);
		}
	}

	for(int j = 0; j < W; j++)
	{
		 g[i * W + j] += log10(pow(re[j], 2.0f) + pow(im[j], 2.0f));
		Re[i * W + j] += re[j];
		Im[i * W + j] += im[j];
	}

	free(re);
	free(im);
}

void dft_1d_h(float* f, float* g, float* Re, float* Im, int j)
{
	float* re = (float*)calloc(H, sizeof(float));
	float* im = (float*)calloc(H, sizeof(float));

	for(int i = 0; i < H; i++)
	{
		int v = i - H / 2;
		for(int y = 0; y < H; y++)
		{
			re[i] +=  f[y * W + j] * cos(2 * M_PI * v * y / H);
			im[i] += -f[y * W + j] * sin(2 * M_PI * v * y / H);
		}
	}

	for(int i = 0; i < H; i++)
	{
		 g[i * W + j] += log10(pow(re[i], 2.0f) + pow(im[i], 2.0f));
		Re[i * W + j] += re[i];
		Im[i * W + j] += im[i];
	}
	free(re);
	free(im);

}


void IDFT(float* Re, float* Im, unsigned char fname2[])
{
	float* g = (float*)calloc(H * W, sizeof(float));

	for(int i = 0; i < H; i++)
	{
		idft_1d_w(Re, Im, g, i);
	}

	for(int j = 0; j < W; j++)
	{
		idft_1d_h(Re, Im, g, j);
	}

	for(int j = 0; j < W*W/2; j++)
	{
		printf("g[%d] = %f\n", j, g[j]);
	}

	writeRawFile(fname2, sizeof(float), H * W, g);

	free(g);
}

void idft_1d_w(float* Re, float* Im, float* g, int i)
{
	float* re = (float*)calloc(W, sizeof(float));


	for(int x = 0; x < W; x++)
	{
		for(int j = 0; j < W; j++)
		{
			int u = j - W / 2;
			re[x] +=  Re[i * W + j] * cos(2 * M_PI * u * x / W) / W - Im[i * W + j] * sin(2 * M_PI * u * x / W) / W;
		}
	}

	for(int j = 0; j < W; j++)
	{
		g[i * W + j]  += re[j];
	}

	free(re);

}

void idft_1d_h(float* Re, float* Im, float* g, int j)
{
	float* re = (float*)calloc(H, sizeof(float));


	for(int y = 0; y < H; y++)
	{
		for(int i = 0; i < H; i++)
		{
			int v = i - H / 2;
			re[y] +=  Re[i * W + j] * cos(2 * M_PI * v * y / H) / H - Im[i * W + j] * sin(2 * M_PI * v * y / H) / H;
		}
	}

	for(int i = 0; i < H; i++)
	{
		g[i * W + j]  += re[i];
	}


	free(re);
}






void readRawFile (const unsigned char fname[], const size_t size, const size_t num, float* image)
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
