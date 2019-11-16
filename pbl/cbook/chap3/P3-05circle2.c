/* P3-05circle2.c */

#include  <stdio.h>
#include  <stdlib.h>
#define  N  4

main( )
{
	char    fi[50];
	float   img[128*128];
	int     nx = 128, ny = 128;
	int     i, j, kx, ky;
	double  x, y, x0, y0, r = 32.;
	FILE    *fp;

	// ‰æ‘œ‚Ì‰Šú‰»
	for (i = 0 ; i < nx*ny ; i++)
		img[i] = 0;

	// ‰~‰æ‘œ‚Ìì¬
	for (i = 0 ; i < ny ; i++) {
		y0 = ny/2-i;
		for(j = 0 ; j < nx ; j++) {
			x0 = j-nx/2;
			for(ky = 0 ; ky < N ; ky++) {
				y = y0+(N-2*ky-1)/(double)(2*N);
				for(kx = 0 ; kx < N ; kx++) {
					x = x0+(2*kx+1-N)/(double)(2*N);
					if(x*x+y*y <= r*r)
						img[i*nx+j] += 100;
				}
			}
			img[i*nx+j] /= N*N;
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
