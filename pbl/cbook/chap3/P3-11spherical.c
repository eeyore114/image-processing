/* P3-11spherical.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  10     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char    f1[50]; /* output new image file name */
	char    f2[50]; /* output new image file name */
	float   *simg;  /* sin image matrix */
	float   *cimg;  /* cos image matrix */
	int     nx;     /* number of width */
	int     ny;     /* number of height */
	double  a;      /* amplitude of wave */
	double  u;      /* frequency */
	double  x0;     /* center of x-direction */
	double  y0;     /* center of y-direction */
	double  p;      /* initial phase */
} Param;

char *menu[PN] = {
	"sine and cosine spherical wave",
	"Output sin image file name <float> ",
	"Output cos image file name <float> ",
	"Number of width       ",
	"Number of height      ",
	"amplitude of wave     ",
	"frequency of wave     ",
	"center of x-direction ",
	"center of y-direction ",
	"initial phase         ",
};

void write_data(char *, float *, int);
void spherical_wave(float *, float *, int, int, double, double, double, double, double);

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
	sprintf( pm->f1, "ssin.img");
	sprintf( pm->f2, "scos.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->a  = 1;
	pm->u  = 5;
	pm->x0 = 0;
	pm->y0 = 0;
	pm->p  = 0;

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
		fprintf( stdout, " %s [%f] :", menu[i++], pm->a );
		if(*gets(dat) != '\0')  pm->a = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->u );
		if(*gets(dat) != '\0')  pm->u = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->x0 );
		if(*gets(dat) != '\0')  pm->x0 = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->y0 );
		if(*gets(dat) != '\0')  pm->y0 = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->p );
		if(*gets(dat) != '\0')  pm->p = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->a  = atof( argv[i++] );
		if((argc--) > 1) pm->u  = atof( argv[i++] );
		if((argc--) > 1) pm->x0 = atof( argv[i++] );
		if((argc--) > 1) pm->y0 = atof( argv[i++] );
		if((argc--) > 1) pm->p  = atof( argv[i++] );
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

	pm->simg = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->cimg = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Make Sine & Cosine spherical wave Image ***\n");
	spherical_wave(pm->simg, pm->cimg, pm->nx, pm->ny, pm->a, pm->u, pm->x0, pm->y0, pm->p);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f1, pm->simg, pm->nx*pm->ny);
	write_data(pm->f2, pm->cimg, pm->nx*pm->ny);

	free(pm->simg);
	free(pm->cimg);
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

void spherical_wave(float *simg, float *cimg, int nx, int ny, double a, double u, double x0, double y0, double p)
{
	int    i, j;
	double x, y;

	for(i = 0 ; i < ny ; i++) {
		y = ny/2-i;
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;
			simg[i*nx+j] = (float)(a*sin(2*PI*u*sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/nx+p));
			cimg[i*nx+j] = (float)(a*cos(2*PI*u*sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/nx+p));
		}
	}
}
