#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#define W 256
#define H 256

const  char readFileName[] = "Brain_uchar_256-256.raw";
const  char writeFileName[] = "Brain_float_256-256.raw";


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image);

void writeRawFile(const char fname[], const size_t size, const size_t num, float* image);

template <typename T> 
T add(T a, T b);

template <typename T> 
T add(T a, T b) { return (a + b); }


int main()
{
	float* g = (float*)calloc(H * W, sizeof(float));
	{
		unsigned char* f = (unsigned char*)calloc(H * W, sizeof(unsigned char));
	/*
		FILE* fp1 = fopen(readFileName , "rb");
		fread( f , sizeof(unsigned char) , W * H , fp1);
		fclose( fp1 );
	*/
		readRawFile(readFileName, H * W, f);	
		for(int i = 0; i < W * H; i++)
		{
			g[i] = (float)f[i];
		}
	}

	writeRawFile(writeFileName, sizeof(float), H * W, g);

	{
		float a = 1.0, b = 2.0;
		std::cout << add(a, b) << std::endl;
	}

	{
		int a = 3, b = 2;
		std::cout << add(a, b) << std::endl;
	}

	{
		std::string a = "abe", b = "konoshin";
		std::cout << add(a, b) << std::endl;
	}

}


template <typename T> void 
readRawFile (const char fname[], const size_t num, T* image)
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
