/*	P4-17spectrum.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  6     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input real      file name */
	char   f2[50]; /* input imaginary file name */
	char   f3[50]; /* output file name */
	float  *img;   /* real image */
	float  *imi;   /* imaginary image */
	float  *spt;   /* spectrum */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
} Param;

char *menu[PN] = {
	"Power spectrum of image",
	"Input real      file name <float> ",
	"Input imaginary file name <float> ",
	"Output file name          <float> ",
	"Number of matrix  (x)    ",
	"Number of matrix  (y)    ",
};

void read_data(char *, float *, int);
void write_data_ascii(char *, float *, int);
void spectrum(float *, float *, float *, int, int);

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
	sprintf( pm->f1, "n1r.img");
	sprintf( pm->f2, "n1i.img");
	sprintf( pm->f3, "s0.txt");
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

	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->imi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->spt = (float *)malloc((unsigned long)pm->nx/2*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->img, pm->nx*pm->ny);
	read_data(pm->f2, pm->imi, pm->nx*pm->ny);

	printf(" *** Calculation spectrum  ***\n");
	spectrum(pm->spt, pm->img, pm->imi, pm->nx, pm->ny);

	printf(" *** Write Spectrum data   ***\n");
	write_data_ascii(pm->f3, pm->spt, pm->nx/2);

	free(pm->img);
	free(pm->imi);
	free(pm->spt);
	free(pm);
}

void read_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and read data */
	if(NULL == (fp = fopen(fi, "rb"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(prj, sizeof(float), size, fp);
	fclose(fp);
}

void write_data_ascii(char *fi, float *prj, int size)
{
	int		i;
	FILE   *fp;

	/* open file and write data */
	if(NULL == (fp = fopen(fi, "w"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	for(i = 0 ; i < size ; i++) {
		fprintf(fp, "%f\n", prj[i]);
	}
	fclose(fp);
}

void spectrum(float *spt, float *img, float *imi, int nx, int ny)
{
	int     i, j, ix, iy, *n;
	double  x, y, t, dx, dy, re, im;

	n = (int *)malloc((unsigned long)nx*sizeof(int));
	for(i = 0 ; i < nx/2 ; i++) {
		spt[i] = 0;
		n[i] = 0;
	}
	for(i = 0 ; i < nx/2 ; i++) {
		for(j = 0 ; j < ny ; j++) {
			t = 2*j*PI/ny;
			x = i*cos(t);
			y = i*sin(t);
			ix = (int)(x+nx/2);
			dx = x+nx/2-ix;
			iy = (int)(ny/2-y);
			dy = ny/2-y-iy;
			if(ix < 0 || ix > nx-2 || iy < 0 || iy > ny-2) continue;
			re = (1-dy)*(1-dx)*img[iy*nx+ix]+(1-dy)*dx*img[iy*nx+ix+1]
			    +dy*(1-dx)*img[(iy+1)*nx+ix]+dy*dx*img[(iy+1)*nx+ix+1];
			im = (1-dy)*(1-dx)*imi[iy*nx+ix]+(1-dy)*dx*imi[iy*nx+ix+1]
			    +dy*(1-dx)*imi[(iy+1)*nx+ix]+dy*dx*imi[(iy+1)*nx+ix+1];
			spt[i] += (float)(re*re+im*im);
			n[i]++;
		}
	}
	for(i = 0 ; i < nx/2 ; i++)
		if(n[i] != 0)
			spt[i] /= n[i];
	free(n);
}
