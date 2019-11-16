
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#include "xcom.h"
#define detecter_size 65.

int size = (int)detecter_size;
int photon_num = 100000000;

template <class T> 
void readRawFile (const char fname[], const size_t num, T* image);

template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);

/*
やること



*/

int main()
{

}


template <class T> 
void readRawFile (const char fname[], const size_t num, T* image)
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




template <class T>
void writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}

