#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#define PI 3.1415
#define M  360

char readfilename[] 	= "datasin.csv";
char writefilename[] 	= "sin_dft.csv";

char dft_1d(char fname1[], char fname2[]);


void main(void)
{
	dft_1d(readfilename, writefilename);
}




char dft_1d(char fname1[], char fname2[])
{
	FILE *fp1, *fp2;
	double* re = (double*)calloc(M, sizeof(double));
	double* im = (double*)calloc(M, sizeof(double));
	double* f  = (double*)calloc(M, sizeof(double));

	fp1 = fopen(fname1, "r");
	if(fp1 == NULL)
	{
		printf("file1 not open");
		return 1;
	}

	fp2 = fopen(fname2, "w");
	if(fp2 == NULL)
	{
		printf("file2 not open");
		return 1;
	}

	
	for(int i = 0; i < M; i++)
		fscanf(fp1, "%lf", &f[i]);

	for(int i = 0; i < M; i++)
		printf("%f\n", f[i]);


	for(int i = 0; i < M; i++)
	{
		int u = i - M / 2;
		for(int x = 0; x < M; x++)
		{
			re[i] +=  f[x] * cos(2 * PI * u * x / M);
			im[i] += -f[x] * sin(2 * PI * u * x / M);
		}
		fprintf(fp2, "%d, %f, %f\n", i, re[i], im[i]);
	}

	fclose(fp1);
	fclose(fp2);
	free(re);
	free(im);
	free(f);
	return 0;

}

