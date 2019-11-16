/* P4-24amphesimg.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  7     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char    f1[50]; /* input real part file name */
	char    f2[50]; /* output real part file name */
	char    f3[50]; /* output imaginary part file name */
	float   *fr;    /* real part data */
	float   *fi;    /* imaginary part data */
	float   *ar;    /* amplitude data (real) */
	float   *ai;    /* amplitude data (imaginary) */
	float   *pr;    /* phase data (real) */
	float   *pi;    /* phase data (imaginary) */
	int     nx;     /* number of width  */
	int     ny;     /* number of height */
	double  th;     /* threshold */
} Param;

char *menu[PN] = {
	"Amplitude and Phase image (Real space)",
	"Input  image     file name <float> ",
	"Output amplitude file name <float> ",
	"Output phase     file name <float> ",
	"Number of width    ",
	"Number of height   ",
	"Threshold          ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void fft2d(int, float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void amplitude(float *, float *, float *, float *, int);
void phase(float *, float *, float *, float *, int, double);

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
	sprintf( pm->f2, "n1a.img");
	sprintf( pm->f3, "n1p.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->th = 0.01;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f3 );
		if(*gets(dat) != '\0')  strcpy(pm->f3, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->th );
		if(*gets(dat) != '\0')  pm->th = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->th = atof( argv[i++] );
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

	pm->fr = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->fi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->ar = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->ai = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->pr = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->pi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fr, pm->nx*pm->ny);
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->fi[i] = 0;

	printf(" *** %s ***\n", menu[0]);
	fft2d(1, pm->fr, pm->fi, pm->nx, pm->ny);

	// U•‰æ‘œ‚ÌŒvŽZ
	amplitude(pm->ar, pm->ai, pm->fr, pm->fi, pm->nx*pm->ny);

	// ˆÊ‘Š‰æ‘œ‚ÌŒvŽZ
	phase(pm->pr, pm->pi, pm->fr, pm->fi, pm->nx*pm->ny, pm->th);

	fft2d(-1, pm->ar, pm->ai, pm->nx, pm->ny);
	fft2d(-1, pm->pr, pm->pi, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->ar, pm->nx*pm->ny);
	write_data(pm->f3, pm->pr, pm->nx*pm->ny);

	free(pm->fr);
	free(pm->fi);
	free(pm->ar);
	free(pm->ai);
	free(pm->pr);
	free(pm->pi);
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

void amplitude(float *ar, float *ai, float *xr, float *xi, int n)
{
	int		i;

	for (i = 0; i < n ; i++) {
		ar[i] = (float)sqrt((double)(xr[i]*xr[i] + xi[i]*xi[i]));
		ai[i] = 0;
	}
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
