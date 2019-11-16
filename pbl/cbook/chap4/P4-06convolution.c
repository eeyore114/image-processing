/* P4-06convolution.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  5     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input  f(x) file name */
	char   f2[50]; /* input  h(x) file name */
	char   f3[50]; /* output g(x) file name */
	float  *fx;    /* f(x) data */
	float  *hx;    /* h(x) data */
	float  *gx;    /* g(x) data */
	int    nx;     /* number of data */
} Param;

char *menu[PN] = {
	"1-D convolution",
	"Input  f(x) file name <txt> ",
	"Input  h(x) file name <txt> ",
	"Output g(x) file name <txt> ",
	"Number of data           ",
};

void read_data_txt(char *, float *, int);
void write_data_txt(char *, float *, int);
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
	sprintf( pm->f1, "fx.txt");
	sprintf( pm->f2, "hx.txt");
	sprintf( pm->f3, "gx.txt");
	pm->nx = 128;

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
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
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

	pm->fx = (float *)malloc((unsigned long)pm->nx*sizeof(float));
	pm->hx = (float *)malloc((unsigned long)pm->nx*sizeof(float));
	pm->gx = (float *)malloc((unsigned long)pm->nx*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data_txt(pm->f1, pm->fx, pm->nx);
	read_data_txt(pm->f2, pm->hx, pm->nx);

	printf(" *** 1-D convolution   ***\n");
	conv1d(pm->gx, pm->fx, pm->hx, pm->nx);

	printf(" *** Write Image data   ***\n");
	write_data_txt(pm->f3, pm->gx, pm->nx);

	free(pm->fx);
	free(pm->hx);
	free(pm->gx);
	free(pm);
}

void read_data_txt(char *fi, float *img, int size)
{
	int    i;
	char   buff[256];
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "r")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	for(i = 0 ; i < size ; i++) {
		fgets(buff, 256, fp);
		img[i] = (float)atof(buff);
	}
	fclose(fp);
}

void write_data_txt(char *fi, float *img, int size)
{
	int    i;
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "w")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	for(i = 0 ; i < size ; i++) {
		fprintf(fp, "%f\n", img[i]);
	}
	fclose(fp);
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
