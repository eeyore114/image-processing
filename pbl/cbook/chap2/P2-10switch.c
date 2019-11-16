/* P2-10switch.c */

#include  <stdio.h>

main( )
{
	char    a;
	printf("Input an alphabet : ");
	scanf("%c", &a);
	switch (a) {
	case   'A':
		printf("This is an apple.\n");
		break;
	case   'B':
		printf("This is a book.\n");
		break;
	case   'C':
		printf("This is a car.\n");
		break;
	default:
		printf("What is this ?\n");
	}
}
