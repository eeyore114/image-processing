#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <random>



template <class T>
void readRawFile (std::string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname.c_str());
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}

template <class T>
void writeRawFile(std::string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname.c_str());
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}
