/* P3-18calc.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  7     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input image file name (1) */
	char   f2[50]; /* input image file name (2) */
	char   f3[50]; /* output new image file name */
	float  *img1;  /* image matrix 1 */
	float  *img2;  /* image matrix 2 */
	int    nx;     /* number of width */
	int    ny;     /* number of height */
	int    tc;     /* type of culculation */
} Param;

char *menu[PN] = {
	"calculation of 2 images",
	"Input  image file name (1) <float> ",
	"Input  image file name (2) <float> ",
	"Output image file name     <float> ",
	"Number of width       ",
	"Number of height      ",
	"Calculation (1:+, 2:-, 3:*, 4:/ )  ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void calculation(float *, float *, int, int, int);

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
	sprintf( pm->f3, "n2.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->tc = 1;

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
		fprintf( stdout, " %s [%d] :", menu[i++], pm->tc );
		if(*gets(dat) != '\0')  pm->tc = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->tc = atoi( argv[i++] );
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

	pm->img1 = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->img2 = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->img1, pm->nx*pm->ny);
	read_data(pm->f2, pm->img2, pm->nx*pm->ny);

	printf(" *** Caluculation of 2 images ***\n");
	calculation(pm->img1, pm->img2, pm->nx, pm->ny, pm->tc);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->img1, pm->nx*pm->ny);

	free(pm->img1);
	free(pm->img2);
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

void calculation(float *img1, float *img2, int nx, int ny, int tc)
{
	int   i;

	switch(tc) {
		case  1:
			for(i = 0 ; i < nx*ny ; i++)
				img1[i] += img2[i];
			break;
		case  2:
			for(i = 0 ; i < nx*ny ; i++)
				img1[i] -= img2[i];
			break;
		case  3:
			for(i = 0 ; i < nx*ny ; i++)
				img1[i] *= img2[i];
			break;
		case 4:
			for(i = 0 ; i < nx*ny ; i++) {
				if(img2[i] != 0.)
					img1[i] /= img2[i];
			}
		}
}
