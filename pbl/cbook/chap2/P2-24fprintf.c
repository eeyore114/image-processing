/* P2-24fprintf.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	FILE   *fp;

	printf( "Input new text file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "w")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fprintf(fp, "Hello C World !! \n");
	fclose (fp);
}
