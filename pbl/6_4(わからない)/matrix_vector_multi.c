#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 2

int main()
{
	int i;
	int j;
	double A[N];
	double B[N*N];
	double C[N];


	for (j = 0; j < N;j++)
	{
		A[j] = 0.0;
	}
	

	for(j = 0;j<N;j++)
	{
		C[j] = 1.0;
	}

	for(i=0;i<N;i++)
	{
		for(j = 0;j<N;j++)
			B[N*j + i] = (double)i/N;
	}
	

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{

			A[i] += B[N*j+i]*C[j];
			 
		}


	}

	for (j = 0; j < N;j++)
	{
		printf("A[%d] = %f\n",j,A[j]);
	}
	


}