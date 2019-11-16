/* P3-04circle1.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	float   img[128*128];
	int     nx = 128, ny = 128, r = 32;
	int     i, j, x, y;
	FILE    *fp;

	// ‰æ‘œ‚Ì‰Šú‰»
	for (i = 0 ; i < nx*ny ; i++)
		img[i] = 0;

	// ‰~‰æ‘œ‚Ìì¬
	for (i = 0 ; i < ny ; i++) {
		y = ny/2-i;  // yÀ•W‚Ö‚Ì•ÏŠ·
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;  // xÀ•W‚Ö‚Ì•ÏŠ·
			if(x*x+y*y <= r*r)
				img[i*nx+j] = 100;
		}
	}

	// ‰æ‘œ‚Ì‘‚«o‚µ
	printf( "Input new file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(img, sizeof(float), nx*ny, fp);
	fclose (fp);
}
