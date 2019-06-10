/* P3-14move.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  7     /* number of parameters + 1 */

typedef struct {
	char    f1[50]; /* input image file name */
	char    f2[50]; /* output new image file name */
	float   *img;   /* image matrix */
	int     nx;     /* number of width */
	int     ny;     /* number of height */
	double  dx;     /* movement x-direction */
	double  dy;     /* movement y-direction */
} Param;

char *menu[PN] = {
	"Move an image",
	"Input  image     file name <float> ",
	"Output new image file name <float> ",
	"Number of width      ",
	"Number of height     ",
	"Movement x-direction ",
	"Movement y-direction ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void move(float *, int, int, double, double);

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
	pm->dx = 0;
	pm->dy = 0;

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
		fprintf( stdout, " %s [%f] :", menu[i++], pm->dx );
		if(*gets(dat) != '\0')  pm->dx = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->dy );
		if(*gets(dat) != '\0')  pm->dy = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->dx = atof( argv[i++] );
		if((argc--) > 1) pm->dy = atof( argv[i++] );
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

	printf(" *** Making Phantom Image ***\n");
	move(pm->img, pm->nx, pm->ny, pm->dx, pm->dy);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

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

void move(float *img, int nx, int ny, double dx, double dy)
{
	int    i, j, i0, i1, j0, j1;
	double x0, y0;
	float  *ima;

	ima = (float *)malloc((unsigned long)nx*ny*sizeof(float));
	for(i = 0 ; i < nx*ny ; i++)
		ima[i] = 0;

	for(i = 0 ; i < ny ; i++) {
		y0 = i+dy; // dyの符号は逆にする（y方向の移動）
		i0 = (int)y0;
		i1 = i0+1;
		if(i0 < 0 || i1 > ny-1) continue;
		for(j = 0 ; j < nx ; j++) {
			x0 = j-dx; // （x方向の移動）
			j0 = (int)x0;
			j1 = j0+1;
			if(j0 < 0 || j1 > nx-1) continue;
			ima[i*nx+j] = (float)((j1-x0)*(i1-y0)*img[i0*nx+j0]
			            + (x0-j0)*(i1-y0)*img[i0*nx+j1]
			            + (j1-x0)*(y0-i0)*img[i1*nx+j0]
			            + (x0-j0)*(y0-i0)*img[i1*nx+j1]);
		}
	}
	for(i = 0 ; i < nx*ny ; i++)
		img[i] = ima[i];
	free(ima);
}
