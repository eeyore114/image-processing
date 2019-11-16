/* P2-22fopen.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	FILE   *fp;

	printf( "Input file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "r")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	printf("File open successfully [%s].\n", fi);
	fclose (fp);
}
