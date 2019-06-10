/* P2-12dowhile.c */

#include  <stdio.h>

main( )
{
	int    a, sum, i;
	printf("Input number (a) : ");
	scanf("%d", &a);
	sum = 0;
	i = 1;
	do {
		sum += i;
		i++;
	} while ( i <= a );
	printf("Total : %d\n", sum);
}
