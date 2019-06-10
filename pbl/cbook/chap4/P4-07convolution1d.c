/* P4-07convolution1d.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  6     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input  f(x,y) file name */
	char   f2[50]; /* input  h(x,y) file name */
	char   f3[50]; /* output g(x,y) file name */
	float  *fx;    /* f(x,y) data */
	float  *hx;    /* h(x,y) data */
	float  *gx;    /* g(x,y) data */
	int    nx;     /* number of width */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"1-D convolution (image)",
	"Input  f(x,y) file name <float> ",
	"Input  h(x,y) file name <float> ",
	"Output g(x,y) file name <float> ",
	"Number of width       ",
	"Number of height      ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void conv1dy(float *, float *, float *, int, int);
void conv1d(float *, float *, float *, int);

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
	sprintf( pm->f2, "h0.img");
	sprintf( pm->f3, "g0.img");
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
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->fx = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->hx = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->gx = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fx, pm->nx*pm->ny);
	read_data(pm->f2, pm->hx, pm->nx*pm->ny);

	printf(" *** %s ***\n", menu[0]);
	conv1dy(pm->gx, pm->fx, pm->hx, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->gx, pm->nx*pm->ny);

	free(pm->fx);
	free(pm->hx);
	free(pm->gx);
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

void conv1dy(float *gx, float *fx, float *hx, int nx, int ny)
{
	int   i;
	
	for(i = 0 ; i < ny ; i++) {
		conv1d(gx+i*nx, fx+i*nx, hx+i*nx, nx);
	}
}

void conv1d(float *gx, float *fx, float *hx, int nx)
{
	int    i, j;
	float  *hx2;

	hx2 = (float *)malloc((unsigned long)2*nx*sizeof(float));

	for(i = 0 ; i < nx/2 ; i++) {
		hx2[i] = hx2[nx+i] = hx[i+nx/2];
		hx2[i+nx/2] = hx2[i+3*nx/2] = hx[i];
	}

	for(i = 0 ; i < nx ; i++) {
		gx[i] = 0;
		for(j = 0 ; j < nx ; j++) {
			gx[i] += fx[j]*hx2[i+nx-j];
		}
	}

	free(hx2);
}
