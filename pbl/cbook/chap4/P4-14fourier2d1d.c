/* P4-14fourier2d1d.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  8     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input real part file name */
	char   f2[50]; /* input imaginary part file name */
	char   f3[50]; /* output real part file name */
	char   f4[50]; /* output imaginary part file name */
	float  *fr;    /* real part data */
	float  *fi;    /* imaginary part data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
	int    ir;     /* forward or inverse */
} Param;

char *menu[PN] = {
	"2-D fourier transform",
	"Input  real      part file name <float> ",
	"Input  imaginary part file name <float> ",
	"Output real      part file name <float> ",
	"Output imaginary part file name <float> ",
	"Number of width          ",
	"Number of height         ",
	"Forward(1) or Inverse(-1)",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void fourier2d1d(int, float *, float *, int, int);

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
	sprintf( pm->f1, "n0r.img");
	sprintf( pm->f2, "0_image");
	sprintf( pm->f3, "n1r.img");
	sprintf( pm->f4, "n1i.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->ir = 1;
	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f3 );
		if(*gets(dat) != '\0')  strcpy(pm->f3, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f4 );
		if(*gets(dat) != '\0')  strcpy(pm->f4, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ir );
		if(*gets(dat) != '\0')  pm->ir = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) strcpy( pm->f4, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->ir = atoi( argv[i++] );
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
	if(strcmp(pm->f2, "0_image"))
		read_data(pm->f2, pm->fi, pm->nx*pm->ny);
	else {
		for(i = 0 ; i < pm->nx*pm->ny ; i++)
			pm->fi[i] = 0;
	}
	printf(" *** 2-D fourier transform   ***\n");
	fourier2d1d(pm->ir, pm->fr, pm->fi, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->fr, pm->nx*pm->ny);
	write_data(pm->f4, pm->fi, pm->nx*pm->ny);

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

void fourier1d(int ir, float *fr, float *fi, int nx)
{
	int     i, j, n = 1;
	float   *gr, *gi;
	double  u, x;

	gr = (float *)malloc((unsigned long)nx*sizeof(float));
	gi = (float *)malloc((unsigned long)nx*sizeof(float));
	for(i = 0 ; i < nx ; i++) {
		u = i-nx/2;
		gr[i] = gi[i] = 0;
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;
			gr[i] += (float)( fr[j]*cos(2*PI*u*x/nx)+ir*fi[j]*sin(2*PI*u*x/nx));
			gi[i] += (float)(-ir*fr[j]*sin(2*PI*u*x/nx)+fi[j]*cos(2*PI*u*x/nx));
		}
	}
	if(ir == -1)  n = nx; // 逆変換はデータ数で割る
	for(i = 0 ; i < nx ; i++) {
		fr[i] = gr[i]/n;
		fi[i] = gi[i]/n;
	}
	free(gr);
	free(gi);
}

void fourier2d1d(int ir, float *fr, float *fi, int nx, int ny)
{
	int   i, j;
	float *gr, *gi;

	// x方向のフーリエ変換
	gr = (float *)malloc((unsigned long)nx*sizeof(float));
	gi = (float *)malloc((unsigned long)nx*sizeof(float));
	for(i = 0 ; i < ny ; i++) {
		fprintf( stderr, "\r   calculation -x- [%3d/%3d]", i+1, ny);
		for(j = 0 ; j < nx ; j++) {
			gr[j] = fr[i*nx+j];
			gi[j] = fi[i*nx+j];
		}
		fourier1d(ir, gr, gi, nx);
		for(j = 0 ; j < nx ; j++) {
			fr[i*nx+j] = gr[j];
			fi[i*nx+j] = gi[j];
		}
	}
	fprintf( stderr, "\n");
	free(gr);
	free(gi);

	// y方向のフーリエ変換
	gr = (float *)malloc((unsigned long)ny*sizeof(float));
	gi = (float *)malloc((unsigned long)ny*sizeof(float));
	for(j = 0 ; j < nx ; j++) {
		fprintf( stderr, "\r   calculation -y- [%3d/%3d]", j+1, nx);
		for(i = 0 ; i < ny ; i++) {
			gr[i] = fr[i*nx+j];
			gi[i] = fi[i*nx+j];
		}
		fourier1d(ir, gr, gi, ny);
		for(i = 0 ; i < ny ; i++) {
			fr[i*nx+j] = gr[i];
			fi[i*nx+j] = gi[i];
		}
	}
	fprintf( stderr, "\n");
	free(gr);
	free(gi);
}
