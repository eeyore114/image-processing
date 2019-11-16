/* P2-21function5.c */

#include  <stdio.h>

void  circle5(double, double *);

main( )
{
	double  a = 10.0, b[2];
	circle5(a, b);
	printf("S=%f, L=%f \n", b[0], b[1]);
}

void  circle5(double  r, double  *s)
{
	s[0] = 3.14*r*r;
	s[1] = 2*3.14*r;
}
