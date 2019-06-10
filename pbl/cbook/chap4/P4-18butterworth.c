/*  P4-18butterworth.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  7

typedef struct {
	char    f1[50]; /* input file name */
	char    f2[50]; /* output file name */
	float   *img;   /* real image */
	float   *imi;   /* imaginary image */
	int     nx;     /* number of matrix (x) */
	int     ny;     /* number of matrix (y) */
	double  ct;     /* cutoff of butterworth */
	int     od;     /* order of butterworth */
} Param;

char *menu[PN] = {
	"Butterworth filter",
	"Input  file name <float> ",
	"Output file name <float> ",
	"Number of matrix  (x)    ",
	"Number of matrix  (y)    ",
	"Cutoff of Butterworth    ",
	"Order  of Butterworth    ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void butterworth(float *, int, int, double, int);
void fft2d(int, float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);

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
	pm->ct = 0.25;
	pm->od = 4;

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
		fprintf( stdout, " %s [%5.2f] :", menu[i++], pm->ct );
		if(*gets(dat) != '\0')  pm->ct = atof(dat);
		fprintf( stdout, " %s [%5d] :", menu[i++], pm->od );
		if(*gets(dat) != '\0')  pm->od = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->ct = atof( argv[i++] );
		if((argc--) > 1) pm->od = atoi( argv[i++] );
	}
	else {
		usage(argc, argv);
	}
}

main(int argc, char *argv[] )
{
	int     i;
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->imi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->img, pm->nx*pm->ny);
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->imi[i] = 0;

	printf(" *** 2-D fourier transform   ***\n");
	fft2d(1, pm->img, pm->imi, pm->nx, pm->ny);

	printf(" *** Butterworth filtering  ***\n");
	butterworth(pm->img, pm->nx, pm->ny, pm->ct, pm->od);
	butterworth(pm->imi, pm->nx, pm->ny, pm->ct, pm->od);

	printf(" *** 2-D inverse fourier transform   ***\n");
	fft2d(-1, pm->img, pm->imi, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm->imi);
	free(pm);
}

void read_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and read data */
	if((fp = fopen(fi, "rb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(prj, sizeof(float), size, fp);
	fclose(fp);
}

void write_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(prj, sizeof(float), size, fp);
	fclose(fp);
}

void butterworth(float *img, int nx, int ny, double ct, int od)
{
	int     i, j;
	float   *flt;
	double  fx, fy, u;

	flt = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	for(i = 0 ; i < ny ; i++) {
		fy = 1.0*(ny/2-i)/ny;
		for(j = 0 ; j < nx ; j++) {
			fx = 1.0*(j-nx/2)/nx;
			u = sqrt(fx*fx+fy*fy);
			flt[i*nx+j] = (float)(1./(1.+pow(u/ct, 2.0*od)));
			img[i*nx+j] *= flt[i*nx+j];
		}
	}
	write_data("filter.img", flt, nx*ny);
	free(flt);
}
