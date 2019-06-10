#include <iostream>
#include <math.h>
using namespace std;

#define W 256
#define H 256

const char readfilename1[] = "Brain_float_512-512.raw";
const char readfilename2[] = "ml-em_100_float_512-512.raw";


void readRawFile (const char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);


int main()
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	float sum_signal = 0.0f;
	float sum_noise  = 0.0f;

	readRawFile(readfilename1, sizeof(float), H * W, f);
	readRawFile(readfilename2, sizeof(float), H * W, g);

	for(int i = 0; i < H; i++)
	{
		for(int j = 0; j < W; j++)
		{
			sum_signal += pow(f[i * W + j], 2.0f);
			sum_noise  += pow(f[i * W + j] - g[i * W + j], 2.0f); 
		}
	}


	float snr = 10 * log10f(sum_signal / sum_noise);

	cout << "SNR = " << snr << "\n";


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