/* P3-06ellipse1.c */

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#define  PI  3.14159265358979
#define  N   4

main( )
{
	char    fi[50];
	float   img[128*128];
	int     nx = 128, ny = 128;
	int     i, j, kx, ky;
	double  x0 = 10., y0 = 10., a = 20., b = 30., th = 15.;
	double  x, y, x1, y1, x2, y2, si, co;
	FILE    *fp;

	// ‰æ‘œ‚Ì‰Šú‰»
	for (i = 0 ; i < nx*ny ; i++)
		img[i] = 0;

	// ‘È‰~‰æ‘œ‚Ìì¬
	th = -th*PI/180.;
	si = sin(th);
	co = cos(th);
	for (i = 0 ; i < ny ; i++) {
		y1 = ny/2-i-y0;
		for(j = 0 ; j < nx ; j++) {
			x1 = j-nx/2-x0;
			for(ky = 0 ; ky < N ; ky++) {
				y2 = y1+(N-2*ky-1)/(double)(2*N);
				for(kx = 0 ; kx < N ; kx++) {
					x2 = x1+(2*kx+1-N)/(double)(2*N);
					x = x2*co-y2*si;
					y = x2*si+y2*co;
					if(x*x/(a*a)+y*y/(b*b) <= 1.)
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
