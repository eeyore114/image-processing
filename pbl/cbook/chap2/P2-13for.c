/* P2-13for.c */

#include  <stdio.h>

main( )
{
	int    a, sum, i;
	printf("Input number (a) : ");
	scanf("%d", &a);
	sum = 0;
	for ( i = 1 ; i <= a ; i++ ) {
		sum += i;
	}
	printf("Total : %d\n", sum);
}
