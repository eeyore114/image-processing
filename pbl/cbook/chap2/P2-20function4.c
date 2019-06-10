/* P2-20function4.c */

#include  <stdio.h>

void circle4(double *);

main( )
{
	double  a = 10.0;
	circle4(&a);
	printf("a=%f \n", a);
}

void  circle4(double *r)
{
	*r = 3.14*(*r)*(*r);
}
