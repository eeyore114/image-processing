/* P4-08convolution2d.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  6     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input  f(x,y) file name */
	char   f2[50]; /* input  h(x,y) file name */
	char   f3[50]; /* output g(x,y) file name */
	float  *fxy;    /* f(x,y) data */
	float  *hxy;    /* h(x,y) data */
	float  *gxy;    /* g(x,y) data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"2-D convolution",
	"Input  f(x,y) file name <float> ",
	"Input  h(x,y) file name <float> ",
	"Output g(x,y) file name <float> ",
	"Number of width   ",
	"Number of height  ",
	};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void conv2d(float *, float *, float *, int, int);

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

	pm->fxy = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->hxy = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->gxy = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fxy, pm->nx*pm->ny);
	read_data(pm->f2, pm->hxy, pm->nx*pm->ny);

	printf(" *** 2-D convolution   ***\n");
	conv2d(pm->gxy, pm->fxy, pm->hxy, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->gxy, pm->nx*pm->ny);

	free(pm->fxy);
	free(pm->hxy);
	free(pm->gxy);
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

void conv2d(float *gxy, float *fxy, float *hxy, int nx, int ny)
{
	int   i, j, m, n;
	float *hxy2, buff;

	// 重畳積分のため、縦横2倍にする
	hxy2 = (float *)malloc((unsigned long)4*nx*ny*sizeof(float));
	// 縦横2倍にする画像を4分割で反転させる
	for(i = 0 ; i < ny/2 ; i++) {
		for(j = 0 ; j < nx/2 ; j++) {
			buff = hxy[i*nx+j];
			hxy[i*nx+j] = hxy[(i+ny/2)*nx+j+nx/2];
			hxy[(i+ny/2)*nx+j+nx/2] = buff;
			buff = hxy[i*nx+j+nx/2];
			hxy[i*nx+j+nx/2] = hxy[(i+ny/2)*nx+j];
			hxy[(i+ny/2)*nx+j] = buff;
		}
	}
	// 反転画像を縦横2倍画像の4カ所に当てはめる
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			hxy2[i*nx*2+j] = hxy2[i*nx*2+j+nx] = hxy2[(i+ny)*nx*2+j] = hxy2[(i+ny)*nx*2+j+nx] = hxy[i*nx+j];
		}
	}

	// 2次元重畳積分の実行
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			gxy[i*nx+j] = 0;
			for(m = 0 ; m < ny ; m++) {
				for(n = 0 ; n < nx ; n++) {
					gxy[i*nx+j] += fxy[m*nx+n]*hxy2[(i+ny-m)*nx*2+j+nx-n];
				}
			}
		}
	}

	free(hxy2);
}
