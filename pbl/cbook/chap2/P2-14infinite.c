/* P2-14infinite.c */

#include  <stdio.h>

main( )
{
	int    a, sum, i;
	printf("Input number (a) : ");
	scanf("%d", &a);
	sum = 0;
	i = 1;
	while ( 1 ) {
		sum += i;
		if ( i >= a )  break;
		i++;
	}
	printf("Total : %d\n", sum);
}
