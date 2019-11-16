#include <stdio.h>
#include <math.h>

void rotateCoordinate(float *x,float *y,float theta);
void printCoodinate(float x, float y);


int main()
{
	float x = 1;
	float y = 1;
	float theta = M_PI / 2;
		
	printCoodinate(x, y);
	
	rotateCoordinate(&x,&y,theta);
		
	printCoodinate(x, y);
	
}



void rotateCoordinate(float *x,float *y,float theta)
{	
	float x0 = *x;
	float y0 = *y;
	
	*x = *x * cosf(theta) - *y * sinf(theta);
	*y = *y * cosf(theta) + *x * sinf(theta);	
}

void printCoodinate(float x, float y)
{
	printf("(x, y) = (%3.1f, %3.1f)", x, y);
}
