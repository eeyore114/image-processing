/* P4-04lsf.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  6     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char    f1[50]; /* output new image file name */
	float   *img;   /* image matrix */
	int     nx;     /* number of width */
	int     ny;     /* number of height */
	double  w;      /* FWHM */
	double  x0;     /* center of x-direction */
} Param;

char *menu[PN] = {
	"line spread function image",
	"Output image file name <float> ",
	"Number of width       ",
	"Number of height      ",
	"FWHM                  ",
	"center of x-direction ",
};

void write_data(char *, float *, int);
void lsf(float *, int, int, double, double);

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
	sprintf( pm->f1, "lsf.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->w  = 32;
	pm->x0 = 0;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->w );
		if(*gets(dat) != '\0')  pm->w = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->x0 );
		if(*gets(dat) != '\0')  pm->x0 = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->w  = atof( argv[i++] );
		if((argc--) > 1) pm->x0 = atof( argv[i++] );
	}
	else {
		usage(argc, argv);
	}

}

main(int argc, char *argv[] )
{
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Make lsf Image ***\n");
	lsf(pm->img, pm->nx, pm->ny, pm->w, pm->x0);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f1, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
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

void lsf(float *img, int nx, int ny, double w, double x0)
{
	int     i, j;
	double  a, x;

	a = sqrt(2.7725887/PI/w/w);
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;
			img[i*nx+j] = (float)(a*exp(-2.7725887*(x-x0)*(x-x0)/(w*w)));
		}
	}
}
