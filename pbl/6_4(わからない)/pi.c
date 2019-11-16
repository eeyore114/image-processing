#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 10000

int main()
{
	int i;
	double seed;//乱数
	double x;
	double y;
	int count_r;
	double function;
	double prob;
	clock_t start,end;
	count_r = 0;
	prob = 0.0;
	srand((unsigned)time(NULL));
	int s = 0;



	//時間計測開始
	start = clock();

	for(i = 0;i<N;i++)
	{
	//2次元平面上の点(x,y)を乱数で発生(0<x<1,0<y<1)
	seed = (double)rand()/RAND_MAX;	
	x = seed;
	seed = (double)rand()/RAND_MAX;
	y = seed;

	//printf("x = %f\n", x);
	//printf("y = %f\n", y);

	if(y <= sqrt(1 - pow(x,2)))
		s++;

	
	}
	//時間計測終了
	end = clock();

	printf("count_r = %d\n",count_r);

	prob = (double)count_r/N;

	printf("time = %f sec\n", (double)(end-start)/CLOCKS_PER_SEC);


	printf("S = %f\n", (float)4*s/N);

	return 0;




}