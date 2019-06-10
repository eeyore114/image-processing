/* P3-21profile.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  7

typedef struct {
	char   f1[50]; /* input  file name */
	char   f2[50]; /* profile name */
	float  *img;   /* image matrix */
	int    nx;     /* x size (width) */
	int    ny;     /* y size (height) */
	int    lx;     /* line x point */
	int    ly;     /* line y point */
} Param;

char *menu[PN] = {
	"Make profile for an image",
	"Input  image   file name ",
	"Output profile file name ",
	"x-Size (width)  ",
	"y-Size (height) ",
	"x-point (line)  ",
	"y-point (line)  ",
};

void read_data(char *, float *, int);
void write_profile(char *, float *, int, int, int, int);

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
	sprintf( pm->f2, "n1.txt");
	pm->nx = 128;
	pm->ny = 128;
	pm->lx = 64;
	pm->ly = 64;

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
		fprintf( stdout, " %s [%d] :", menu[i++], pm->lx );
		if(*gets(dat) != '\0')  pm->lx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ly );
		if(*gets(dat) != '\0')  pm->ly = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->lx = atoi( argv[i++] );
		if((argc--) > 1) pm->ly = atoi( argv[i++] );
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

	printf(" *** Write profile data   ***\n");
	write_profile(pm->f2, pm->img, pm->nx, pm->ny, pm->lx, pm->ly);

	free(pm->img);
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

void write_profile(char *fi, float *img, int nx, int ny, int lx, int ly)
{
	int    i;
	FILE   *fp;

	/* open file and write data */
	if(NULL == (fp = fopen(fi, "w"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fprintf(fp,"%s\n", fi);
	fprintf(fp,"x(y=%d)\n", ly);
	for(i = 0 ; i < nx ; i++) {
		fprintf(fp,"%f\n", img[ly*nx+i]);
	}
	fprintf(fp,"\ny(x=%d)\n", lx);
	for(i = 0 ; i < ny ; i++) {
		fprintf(fp,"%f\n", img[i*nx+lx]);
	}
	fclose(fp);
}
