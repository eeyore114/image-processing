/* P2-23fread.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	short   buff[128*128];
	FILE   *fp;

	printf( "Input image file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "rb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fread(buff, sizeof(short), 128*128, fp);
	fclose (fp);

	printf( "Input new file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(buff, sizeof(short), 128*128, fp);
	fclose (fp);
}
