/* P4-21correlation_f.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  6     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input  f(x,y) file name */
	char   f2[50]; /* input  g(x,y) file name */
	char   f3[50]; /* output r(x,y) file name */
	float  *fre;   /* f(x,y) real data */
	float  *fim;   /* f(x,y) imaginary data */
	float  *gre;   /* g(x,y) real data */
	float  *gim;   /* g(x,y) imaginary data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"2-D correlation (Fourier space)",
	"Input  f(x,y) file name <float> ",
	"Input  g(x,y) file name <float> ",
	"Output r(x,y) file name <float> ",
	"Number of width          ",
	"Number of height         ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void fft2d(int, float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void correlation_f(float *, float *, float *, float *, int);

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
	sprintf( pm->f1, "f0.img");
	sprintf( pm->f2, "g0.img");
	sprintf( pm->f3, "r0.img");
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

	pm->fre = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->fim = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->gre = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->gim = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fre, pm->nx*pm->ny);
	read_data(pm->f2, pm->gre, pm->nx*pm->ny);
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->fim[i] = pm->gim[i] = 0;

	printf(" *** %s ***\n", menu[0]);
	fft2d(1, pm->fre, pm->fim, pm->nx, pm->ny);
	fft2d(1, pm->gre, pm->gim, pm->nx, pm->ny);
	correlation_f(pm->fre, pm->fim, pm->gre, pm->gim, pm->nx*pm->ny);
	fft2d(-1, pm->fre, pm->fim, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->fre, pm->nx*pm->ny);

	free(pm->fre);
	free(pm->fim);
	free(pm->gre);
	free(pm->gim);
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

void correlation_f(float *fre, float *fim, float *gre, float *gim, int n)
{
	int   i;
	float rre, rim;

	for(i = 0 ; i < n ; i++) {
		rre = fre[i]*gre[i]+fim[i]*gim[i];
		rim = fim[i]*gre[i]-fre[i]*gim[i];
		fre[i] = rre;
		fim[i] = rim;
	}
}
