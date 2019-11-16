/* P2-09if2.c */

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
	else  if (a > b) {
		printf("a is greater than b.\n");
	}
	else {
		printf("a is less than b.\n");
	}
}
