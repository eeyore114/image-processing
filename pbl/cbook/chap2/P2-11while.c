/* P2-11while.c */

#include  <stdio.h>

main( )
{
	int    a, sum, i;
	printf("Input number (a) : ");
	scanf("%d", &a);
	sum = 0;
	i = 1;
	while ( i <= a ) {
		sum += i;
		i++;
	}
	printf("Total : %d\n", sum);
}
