/* P2-08if1.c */

#include  <stdio.h>

main( )
{
	int	a, b;
	printf("Input number (a) : ");
	scanf("%d", &a);
	printf("Input number (b) : ");
	scanf("%d", &b);
	if (a == b) {
		printf("a is equal to b.\n");
	}
	else {
		printf("a is not equal to b.\n");
	}
}
