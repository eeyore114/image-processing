/* P4-23amphese.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  6     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input real part file name */
	char   f2[50]; /* output real part file name */
	char   f3[50]; /* output imaginary part file name */
	float  *fr;    /* real part data */
	float  *fi;    /* imaginary part data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"Amplitude and Phase image",
	"Input  image     file name <float> ",
	"Output amplitude file name <float> ",
	"Output phase     file name <float> ",
	"Number of width          ",
	"Number of height         ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void fft2d(int, float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void amphase(float *, float *, int);

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
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
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

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fr, pm->nx*pm->ny);
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->fi[i] = 0;

	printf(" *** %s ***\n", menu[0]);
	fft2d(1, pm->fr, pm->fi, pm->nx, pm->ny);
	amphase(pm->fr, pm->fi, pm->nx*pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->fr, pm->nx*pm->ny);
	write_data(pm->f3, pm->fi, pm->nx*pm->ny);

	free(pm->fr);
	free(pm->fi);
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

void amphase(float *xr, float *xi, int n)
{
	int     i;
	double  a;

	for (i = 0; i < n ; i++) {
		a = sqrt((double)(xr[i]*xr[i] + xi[i]*xi[i]));
		xi[i] = (float)atan2((double)xr[i], (double)xi[i]);
		xr[i] = (float)a;
	}
}

