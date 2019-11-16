#include<stdio.h>
#include<stdlib.h>

void main(void)
{
	int *a;

	a = (int*)malloc(sizeof(int) * 256);
	if(a == NULL)
	{
		printf("failed to allocate memory of a.\n");
		exit(1);
	}

	free(a);
} 