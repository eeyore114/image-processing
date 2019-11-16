/* P3-07ellipse2.c */

#include  <stdio.h>
#include  <math.h>
#define  PI  3.14159265358979
#define  N   4

void mkellipse_phantom(float *img, int nx, int ny, double x0, double y0, double a, double b, double th, double de)
{
	int     i, j, kx, ky;
	double  x, y, x1, y1, x2, y2, si, co, value;

	// ‘È‰~‰æ‘œ‚Ì‰ÁŽZ
	th = -th*PI/180.;
	si = sin(th);
	co = cos(th);
	for (i = 0 ; i < ny ; i++) {
		y1 = ny/2-i-y0;
		for(j = 0 ; j < nx ; j++) {
			x1 = j-nx/2-x0;
			value = 0;
			for(ky = 0 ; ky < N ; ky++) {
				y2 = y1+(N-2*ky-1)/(double)(2*N);
				for(kx = 0 ; kx < N ; kx++) {
					x2 = x1+(2*kx+1-N)/(double)(2*N);
					x = x2*co-y2*si;
					y = x2*si+y2*co;
					if(x*x/(a*a)+y*y/(b*b) <= 1.)
						value += de;
				}
			}
			img[i*nx+j] += (float)(value/(N*N));
		}
	}
}
