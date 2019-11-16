#include <stdio.h>

void watosa(int x,int y,int *sum,int *diff);



int main()
{
	int x = 10;
	int y = 5;
	int sum;
	int diff;
	watosa(x, y, &sum, &diff);
	
		
	printf("wa = %d\n", sum);
	printf("sa = %d\n", diff);
	
	
}


void watosa(int x,int y,int *sum,int *diff)
{
	*sum  = x + y;
	*diff = x - y;	
}



