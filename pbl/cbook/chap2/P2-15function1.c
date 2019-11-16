/* P2-15function1.c */

#include  <stdio.h>

void  circle( double  r )
{
	double  s;
	s = 3.14*r*r;
	printf("s=%f \n", s);
}

main( )
{
	double  a = 10.0;
	circle( a );
}
