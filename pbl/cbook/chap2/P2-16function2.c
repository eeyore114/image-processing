/* P2-16function2.c */

#include  <stdio.h>

double  circle2( double  r )
{
	double  s;
	s = 3.14*r*r;
	return  s ;
}

main( )
{
	double  a = 10.0, b;
	b = circle2( a );
	printf("b=%f \n", b);
}
