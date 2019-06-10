#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//のこぎり波の生成(ファイル名, データ数, 最大値, 最小値, データの間隔)
int createSaw(char filename[], int max, double xmax, double xmin, double dx)
{
	FILE *fp;

	fp = fopen(filename, "w");
	if(fp == NULL)
	{
		printf("file not open.\n");
		exit(1);
	}

	{
		int n = 0;
		double y[100];
		for(double x = xmin; x <= xmax; x += dx)
		{
			y[n] = x;
			n++;
		}

		for(int i = 0;i < max; i++)
		{
			for(int m = 0; m < n; m++)
				fprintf(fp, "%f\n",y[m]);
		}
	}

	fclose(fp);
	return 0;
}

int main(void)
{
	createSaw("data1.csv", 5, 0.0, 0.0, 0.1);
}