/* P2-19cast.c */

#include  <stdio.h>

main( )
{
	int    a = 0x12345678;
	int    *p;
	short  *ps;
	p = &a;
	ps = (short *)p;
	printf("%x , %x \n", *ps, *(ps+1));
}
