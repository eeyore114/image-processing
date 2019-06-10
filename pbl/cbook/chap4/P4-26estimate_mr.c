/* P4-26estimate_mr.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  7     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input real part file name */
	char   f2[50]; /* output real part file name */
	char   f3[50]; /* output imaginary part file name */
	float  *f0r;   /* original    image data */
	float  *f1r;   /* move&rotate image data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
	int    dp;     /* a decimal place (rotation) */
	double  th;     /* threshold */
} Param;

char *menu[PN] = {
	"Estimate movement and rotation",
	"Input original      image file name <float> ",
	"Input move & rotate image file name <float> ",
	"Number of width            ",
	"Number of height           ",
	"A decimal place (rotation) ",
	"Threshold                  ",
};

void  read_data(char *, float *, int);
void  write_data(char *, float *, int);
void  fft2d(int, float *, float *, int, int);
void  FFTInit(int, float *, float *, unsigned short *);
void  FFT(int, int, float *, float *, float *, float *, unsigned short *);
void  phase(float *, float *, float *, float *, int, double);
void  correlation_f(float *, float *,float *, float *, float *, float *, int);
void  rotate(float *, float *, int, int, double);
double calc_ang(float *, float *, int, int, int, double);
void  calc_move(float *, float *, int, int, double, double);

void usage(int argc, char **argv)
{
	int   i;

	fprintf( stderr,"\nUSAGE:\n");
	fprintf( stderr,"\nNAME\n");
	fprintf( stderr,"\n  %s - %s\n", argv[0], menu[0]);
	fprintf( stderr,"\nSYNOPSIS\n");
	fprintf( stderr,"\n  %s [-h] parameters...\n", argv[0]);
	fprintf( stderr,"\nPARAMETERS\n");
	for(i = 1 ; i < PN ; i++)
	  fprintf( stderr,"\n %3d. %s\n", i, menu[i]);
	fprintf( stderr,"\n");
	fprintf( stderr,"\nFLAGS\n");
	fprintf( stderr,"\n  -h  Print Usage (this comment).\n");
	fprintf( stderr,"\n");
	exit(1);
}

void getparameter(int argc, char **argv, Param *pm)
{
	int   i;
	char  dat[256];

	/* default parameter value */
	sprintf( pm->f1, "n0.img");
	sprintf( pm->f2, "n1.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->dp = 1;
	pm->th = 0.01;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->dp );
		if(*gets(dat) != '\0')  pm->dp = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->th );
		if(*gets(dat) != '\0')  pm->th = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->dp = atoi( argv[i++] );
		if((argc--) > 1) pm->th = atof( argv[i++] );
	}
	else {
		usage(argc, argv);
	}

}

main(int argc, char *argv[] )
{
	double  ang;
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->f0r = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->f1r = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->f0r, pm->nx*pm->ny);
	read_data(pm->f2, pm->f1r, pm->nx*pm->ny);

	printf(" *** %s ***\n", menu[0]);
	ang = calc_ang(pm->f0r, pm->f1r, pm->nx, pm->ny, pm->dp, pm->th);

	calc_move(pm->f0r, pm->f1r, pm->nx, pm->ny, pm->th, ang);

	free(pm->f0r);
	free(pm->f1r);
	free(pm);
}

void read_data(char *fi, float *img, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "rb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(img, sizeof(float), size, fp);
	fclose(fp);
}

void write_data(char *fi, float *img, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(img, sizeof(float), size, fp);
	fclose(fp);
}

void phase(float *pr, float *pi, float *xr, float *xi, int n, double th)
{
	int     i;
	double  a;

	for (i = 0; i < n ; i++) {
		a = sqrt((double)(xr[i]*xr[i] + xi[i]*xi[i]));
		if(a < th) {
			pr[i] = pi[i] = 0;
		}
		else {
			pr[i] = (float)(xr[i]/a);
			pi[i] = (float)(xi[i]/a);
		}
	}
}

void correlation_f(float *rre, float *rim, float *fre, float *fim, float *gre, float *gim, int n)
{
	int   i;

	for(i = 0 ; i < n ; i++) {
		rre[i] = fre[i]*gre[i]+fim[i]*gim[i];
		rim[i] = fim[i]*gre[i]-fre[i]*gim[i];
	}
}

void rotate(float *rot, float *img, int nx, int ny, double th)
{
	int     i, j, i0, i1, j0, j1;
	double  x0, y0, si, co;

	for(i = 0 ; i < nx*ny ; i++)
		rot[i] = 0;

	si = sin(th*PI/180.);
	co = cos(th*PI/180.);

	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			x0 = (j-nx/2)*co-(ny/2-i)*si+nx/2; // âÒì]Åixï˚å¸Åj
			j0 = (int)x0;
			j1 = j0+1;
			if(j0 < 0 || j1 > nx-1) continue;
			y0 = ny/2-(j-nx/2)*si-(ny/2-i)*co; // âÒì]Åiyï˚å¸Åj
			i0 = (int)y0;
			i1 = i0+1;
			if(i0 < 0 || i1 > ny-1) continue;
			rot[i*nx+j] = (float)((j1-x0)*(i1-y0)*img[i0*nx+j0]
			            + (x0-j0)*(i1-y0)*img[i0*nx+j1]
			            + (j1-x0)*(y0-i0)*img[i1*nx+j0]
			            + (x0-j0)*(y0-i0)*img[i1*nx+j1]);
		}
	}
}

float max_value(float *img, int nx, int ny, int *x, int *y)
{
	int    i, j;
	float  max;

	max = img[0];
	*x = *y = 0;
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			if(img[i*nx+j] > max) {
				max = img[i*nx+j];
				*x = j;
				*y = i;
			}
		}
	}
	return  max;
}

double max_ang(float *f0r, float *p1r, float *p1i, int nx, int ny, double t, double ang, double th, float *max)
{
	int    i, x, y;
	float  m;
	float  *f1r, *f1i, *p0r, *p0i;

	f1r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	f1i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0i = (float *)malloc((unsigned long)nx*ny*sizeof(float));

	rotate(f1r, f0r, nx, ny, t);
	for(i = 0 ; i < nx*ny ; i++)
		f1i[i] = 0;
	fft2d(1, f1r, f1i, nx, ny);
	phase(p0r, p0i, f1r, f1i, nx*ny, th);
	correlation_f(f1r, f1i, p1r, p1i, p0r, p0i, nx*ny);
	fft2d(-1, f1r, f1i, nx, ny);
	m = max_value(f1r, nx, ny, &x, &y);
	if(m > *max) {
		*max = m;
		ang = t;
	}
	free(f1r);
	free(f1i);
	free(p0r);
	free(p0i);
	return  ang;
}

double calc_ang(float *f0r, float *f1r, int nx, int ny, int dp, double th)
{
	int     i, k, x, y;
	float   max;
	float   *f1i, *f2r, *f2i, *p0r, *p0i, *p1r, *p1i;
	double  t, ang, st, dt;

	f1i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	f2r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	f2i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p1r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p1i = (float *)malloc((unsigned long)nx*ny*sizeof(float));

	for(i = 0 ; i < nx*ny ; i++) {
		f2r[i] = f0r[i];
		f2i[i] = f1i[i] = 0;
	}

	fft2d(1, f2r, f2i, nx, ny);
	fft2d(1, f1r, f1i, nx, ny);
	// à ëäâÊëúÇÃåvéZ
	phase(p0r, p0i, f2r, f2i, nx*ny, th);
	phase(p1r, p1i, f1r, f1i, nx*ny, th);

	correlation_f(f2r, f2i, p1r, p1i, p0r, p0i, nx*ny);
	fft2d(-1, f2r, f2i, nx, ny);
	max = max_value(f2r, nx, ny, &x, &y);
	ang = 0;
	free(f2r);
	free(f2i);
	printf(" Ang = %f, max=%f\n", ang, max);

	for(t = 1 ; t < 360. ; t+=10) {
		ang = max_ang(f0r, p1r, p1i, nx, ny, t, ang, th, &max);
	}
	printf(" Ang = %f, max=%f\n", ang, max);

	dt = 10;
	for(k = -1 ; k < dp ; k++) {
		st = ang;
		dt *= 0.1;
		for(t = st-dt*10 ; t < st+dt*10 ; t+=dt) {
			ang = max_ang(f0r, p1r, p1i, nx, ny, t, ang, th, &max);
		}
	printf(" Ang = %f, max=%f\n", ang, max);
	}

	if(ang > 180.) ang -= 360;
	printf("\nAngle    [ %f degree ] \n", ang);

	fft2d(-1, f1r, f1i, nx, ny);

	free(f1i);
	free(p0r);
	free(p0i);
	free(p1r);
	free(p1i);

	return  ang;
}

void calc_move(float *f0r, float *f1r, int nx, int ny, double th, double ang)
{
	int    i, x, y;
	float  max;
	float  *f1i, *f2r, *f2i, *p0r, *p0i, *p1r, *p1i;

	f1i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	f2r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	f2i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p0i = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p1r = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	p1i = (float *)malloc((unsigned long)nx*ny*sizeof(float));

	for(i = 0 ; i < nx*ny ; i++) {
		f2i[i] = f1i[i] = 0;
	}

	fft2d(1, f1r, f1i, nx, ny);
	phase(p1r, p1i, f1r, f1i, nx*ny, th);

	rotate(f2r, f0r, nx, ny, ang);

	fft2d(1, f2r, f2i, nx, ny);
	phase(p0r, p0i, f2r, f2i, nx*ny, th);

	correlation_f(f2r, f2i, p1r, p1i, p0r, p0i, nx*ny);

	fft2d(-1, f2r, f2i, nx, ny);
	max = max_value(f2r, nx, ny, &x, &y);

	printf("Movement [ %d, %d ] (%d, %d)\n", x-nx/2, ny/2-y, x, y);

	free(f1i);
	free(f2r);
	free(f2i);
	free(p0r);
	free(p0i);
	free(p1r);
	free(p1i);
}
