#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define M 360
#define PI 3.1415

char filename[] = "datasin.csv";

int createSin(char arg[]);

int main(void)
{
	createSin(filename);
}

int createSin(char arg[])
{
	FILE *fp;
	double f[M] = {};

	fp = fopen(arg, "w");
	if(fp == NULL)
	{
		printf("failed to open.\n");
		exit(1);
	}	

	for(int i = 0; i < M; i++)
	{
		f[i] = sin(i * PI / 180);
	}

	for(int i = 0; i < M; i++)
	{
		fprintf(fp, "%f\n", f[i]);
	}

	fclose(fp);
	return 0;
}
