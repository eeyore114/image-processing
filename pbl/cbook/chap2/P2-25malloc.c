/* P2-25malloc.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	short  *buff;
	int     nx, ny;
	FILE   *fp;

	printf( "Input image file name: " );
	scanf( "%s", fi );
	printf( "Input image width  : ");
	scanf( "%d", &nx);
	printf( "Input image height : ");
	scanf( "%d", &ny);

	buff = (short *)malloc(nx*ny*sizeof(short));

	if ((fp = fopen ( fi, "rb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fread(buff, sizeof(short), nx*ny, fp);
	fclose (fp);

	printf( "Input new file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(buff, sizeof(short), nx*ny, fp);
	fclose (fp);
	free(buff);
}
