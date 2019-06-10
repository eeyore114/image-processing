/* P2-26math.c */

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#define   PI  3.14159265358979

main( )
{
	char    buff[256];
	double  a, b;

	printf("Input degree : ");
	gets(buff);

	a = atof( buff );
	b = sin( a*PI/180.);
	printf("sin(%f) = %f\n", a, b);
}
