/* P2-05array.c */

#include  <stdio.h>

main( )
{
	int  a[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	int  b;
	printf("Input month :");
	scanf("%d",&b);
	printf("It has %d days.\n", a[b-1]);
}
