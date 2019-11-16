#include <stdio.h>

int addition(int x,int y);
int subtraction(int x,int y);


int addition(int x,int y)
{	
	return x+y;
}

int subtraction(int x,int y)
{
	return x-y;
}



int main()
{
	int x = 10;
	int y = 5;
	
	int wa = addition(x, y);
	int sa = subtraction(x, y);
	
	printf("wa = %d\n", wa);
	printf("sa = %d\n", sa);
 
}
