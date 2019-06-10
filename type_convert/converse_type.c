#include <stdio.h>
#include <stdlib.h>

#define W 128
#define H 128

const  char readFileName[] = "Brain_uchar_128-128.raw";
const  char writeFileName[] = "Brain_float_128-128.raw";


void readRawFile (const char fname[], const size_t size, const size_t num, float* image);
void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);


int main()
{
	float* g = (float*)calloc(H * W, sizeof(float));
	{
		unsigned char* f = (unsigned char*)calloc(H * W, sizeof(unsigned char));
		FILE* fp1 = fopen(readFileName , "rb");
		fread( f , sizeof(unsigned char) , W * H , fp1);
		fclose( fp1 );
		
		for(int i = 0; i < W * H; i++)
		{
			g[i] = (float)f[i];
		}
	}

	
	

	writeRawFile(writeFileName, sizeof(float), H * W, g);
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
