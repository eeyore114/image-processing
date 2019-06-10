#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N 8
#define W 8
#define H 8

unsigned char readfilename[] 	= "lenna_float_256-256.raw";
unsigned char writefilename[] 	= "DCT_8-8.raw";

void dct(float* f, float* g);
float C(int arg);
void zigzag_scan(float* f, float* g);
void quantization_BrightnessSignal(float* f, float* g);
void quantization_ColorDifferenceSignal(float* f, float* g);
void readRawFile (const unsigned char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const unsigned char fname[], const size_t size, const size_t num, float* image);



void main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* l = (float*)calloc(H * W, sizeof(float));

	readRawFile(readfilename, sizeof(float), H * W, f);
	
	
	dct(f, l);

	writeRawFile(writefilename, sizeof(float), W * H, l);

	free(f);
	free(l);
}


void dct(float* f, float* i)
{
	float* g0 = (float*)calloc(N * N, sizeof(float));
	float* g1 = (float*)calloc(N * N, sizeof(float));
	float* g2 = (float*)calloc(N * N, sizeof(float));

	
	for(int u = 0; u < N; u++)
	{
		for(int v = 0; v < N; v++)
		{
			for(int j0 = 0; j0 < N ; j0++)
			{
				for(int k0 = 0; k0 < N ; k0++)
				{
					g0[u * N + v] += C(u) * C(v) / 4 * f[j0 * N + k0] * cos(((2 * j0 + 1) * u * M_PI) / 16) * cos(((2 * k0 + 1) * v * M_PI) / 16);
				}
			}
		}
	}
	quantization_BrightnessSignal(g0, g1);
	zigzag_scan(g1, g2);

	for(int u = 0; u < N; u++)
	{
		for(int v = 0; v < N; v++)
		{
			i[u * N + v] = g0[u * N + v];
			printf("g0[%d] = %f\n", u * N + v,g0[u * N + v]);
			g0[u * N + v] = 0;
		}
	}

	free(g0);
	free(g1);
	free(g2);

}

float C(int arg)
{
	if(arg == 0)
		return 1 / sqrt(2);
	else
		return 1;
}

void zigzag_scan(float* f, float* g)
{
	g[0] = f[0];
	g[1] = f[0 * W + 1];
	g[2] = f[1 * W + 0];
	{			
		int count = 3;
		int i = 2;
		int j = 0;
		int i_max = 2;
		int j_max = 3;
		int i_up = 1;
		int j_up = 1;
		while(count != 36)
		{
			g[count] = f[i * W + j];
			//printf("count = %d	j = %d	i = %d\n", count, j, i);
			if(i_up)
			{
				if(i == i_max)
				{
					i--;
					i_max += 2;
					i_up = 0;
				}
				else
					i++;
			}
			else
			{
				if(i == 0)
					i_up = 1;
				else
					i--;
			}

			if(j_up)
			{
				if(j_max == j)
				{
					j--;
					j_max += 2;
					j_up = 0;
				}
				else
					j++;
			}
			else
			{
				if(j == 0)
					j_up = 1;
				else
					j--;
			}

			count++;
		}

		int i_min = 1;
		int j_min = 2;
		j = 1;
		i = 7;
		g[36] = f[7 * W + 1];
		//printf("count = %d	j = %d	i = %d\n", count, j, i);
		count++;
		i--;
		j++;
		i_up = 0;

		while(count != N * N)
		{
			i_max = 7;
			j_max = 7;
			g[count] = f[i * W + j];
			//printf("count = %d	j = %d	i = %d\n", count, j, i);
			if(i_up)
			{
				if(i == i_max)
				{
					i_up = 0;
				}
				else
					i++;
			}
			else
			{
				if(i == i_min)
				{
					i_min += 2;
					i++;
					i_up = 1;
				}
				else
					i--;
			}

			if(j_up)
			{
				if(j_max == j)
				{
					j_up = 0;
				}
				else
					j++;
			}
			else
			{
				if(j == j_min)
				{
					j_min += 2;
					j++;
					j_up = 1;
				}
				else
					j--;
			}

			count++;

		}
	}
}

void quantization_BrightnessSignal(float* f, float* g)
{

	int table[N * N] = { 16, 11, 10, 16, 24,   40, 51,  61,
		  			     12, 12, 14, 19, 26,   58, 60,  66,
			    		 14, 13, 16, 24, 40,   57, 69,  57,
			    		 14, 17, 22, 29, 51,   87, 80,  62,
  			    		 18, 22, 37, 56, 68,  109, 103, 77,
			  		     24, 36, 55, 64, 81,  104, 113, 92,
			    		 49, 64, 78, 87, 103, 121, 120, 101,
			    		 72, 92, 95, 98, 112, 100, 103, 99};

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			g[i * N + j] = f[i * W + j] / table[i * N + j];

}


void quantization_ColorDifferenceSignal(float* f, float* g)
{
	int table[N * N];

	for(int i = 0; i < N * N; i++)
		table[i] = 99;

	table[0]  = 17;
	table[1]  = 18;
	table[2]  = 24;
	table[3]  = 47;
	table[8]  = 18;
	table[9]  = 21;
	table[10] = 26;
	table[11] = 66;
	table[16] = 24;
	table[17] = 26;
	table[18] = 56;
	table[19] = 66;
	table[24] = 47;
	table[25] = 66;


	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			g[i * N + j] = f[i * W + j] / table[i * N + j];

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

