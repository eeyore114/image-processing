/*  P4-11median.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  5     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input file name */
	char   f2[50]; /* output file name */
	float  *img;   /* image data */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
} Param;

char *menu[PN] = {
	"Median filter",
	"Input  file name <float> ",
	"Output file name <float> ",
	"Number of matrix  (x)    ",
	"Number of matrix  (y)    ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void median_filter(float *, int, int);

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
	sprintf( pm->f2, "n1.img");
	pm->nx = 128;
	pm->ny = 128;

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
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
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

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->img, pm->nx*pm->ny);

	printf(" *** Calculation Median filter  ***\n");
	median_filter(pm->img, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
}

int compare(const void *e1, const void *e2)
{
	if(*(float *)e1 < *(float *)e2)	return -1;
	if(*(float *)e1 > *(float *)e2) return 1;
	return 0;
}

void median_filter(float *img, int nx, int ny)
{
	int     i, j, x, y;
	float   *sub, fil[9];

	sub = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	for(i = 0 ; i < nx*ny ; i++)
		sub[i] = 0;
	for(i = 1 ; i < ny-1 ; i++) {
		for(j = 1 ; j < nx-1 ; j++) {
			for(y = 0 ; y < 3 ; y++)
				for(x = 0 ; x < 3 ; x++)
					fil[y*3+x] = img[(i+y-1)*nx+j+x-1];
			qsort(fil, 9, 4, compare);
			sub[i*nx+j] = fil[5];
		}
	}
	for(i = 0 ; i < ny ; i++)
		for(j = 0 ; j < nx ; j++)
			img[i*nx+j] = sub[i*nx+j];

	free(sub);
}

void read_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and read data */
	if((fp = fopen(fi, "rb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(prj, sizeof(float), size, fp);
	fclose(fp);
}

void write_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(prj, sizeof(float), size, fp);
	fclose(fp);
}
