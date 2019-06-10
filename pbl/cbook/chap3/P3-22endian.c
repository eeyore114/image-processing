/* P3-22endian.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  5     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input image file name */
	char   f2[50]; /* output new image file name */
	short  *img;   /* image matrix */
	int    nx;     /* number of width */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"Change endian",
	"Input  image     file name <short> ",
	"Output new image file name <short> ",
	"Number of width      ",
	"Number of height     ",
	};

void read_data_s(char *, short *, int);
void write_data_s(char *, short *, int);
void endian(short *, int);

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

	pm->img = (short *)malloc((unsigned long)pm->nx*pm->ny*sizeof(short));

	printf(" *** Read Image data   ***\n");
	read_data_s(pm->f1, pm->img, pm->nx*pm->ny);

	printf(" *** %s ***\n", menu[0]);
	endian(pm->img, pm->nx*pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data_s(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
}

void read_data_s(char *fi, short *img, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "rb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(img, sizeof(short), size, fp);
	fclose(fp);
}

void write_data_s(char *fi, short *img, int size)
{
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(img, sizeof(short), size, fp);
	fclose(fp);
}

void endian(short *img, int n)
// Endianの変換
{
	int   i;
	char  e;

	for(i = 0 ; i < n ; i++) { // バイトの交換（スワップ）
		e = *((char *)(img+i));
		*((char *)(img+i)) = *((char *)(img+i)+1);
		*((char *)(img+i)+1) = e;
	}
}

