/* P2-17function3.c */

#include  <stdio.h>

double  r, s;

void  circle3( )
{
	s = 3.14*r*r;
}

main( )
{
	r = 10.0;
	circle3( );
	printf("s=%f \n", s);
}
