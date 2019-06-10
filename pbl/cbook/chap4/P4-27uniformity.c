/* P4-27uniformity.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  4     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input real part file name */
	float  *img;   /* image data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"Integral uniformity and differential uniformity",
	"Input original      image file name <float> ",
	"Number of width            ",
	"Number of height           ",
};

void read_data(char *, float *, int);
void integral_uniformity(float *, int, int);
void differential_uniformity(float *, int, int);

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
	pm->nx = 128;
	pm->ny = 128;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
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

	printf(" *** %s ***\n", menu[0]);
	integral_uniformity(pm->img, pm->nx, pm->ny);
	differential_uniformity(pm->img, pm->nx, pm->ny);

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

void integral_uniformity(float *img, int nx, int ny)
{
	int    i;
	float  max, min, uni;

	// ç≈ëÂílÇ∆ç≈è¨ílÇÃéZèo
	max = min = img[0];
	for(i = 1 ; i < nx*ny ; i++) {
		if(img[i] > max)      max = img[i];
		else if(img[i] < min) min = img[i];
	}
	// ãœàÍê´ÇÃåvéZ
	uni = (max-min)/(max+min)*100;
	printf("\n Integral uniformity = %f\n", uni);
	printf("   [max = %f, min = %f]\n", max, min);
}

void differential_uniformity(float *img, int nx, int ny)
{
	int    i, j, k;
	float  max, min, uni;
	float *dif;

	dif = (float *)malloc((unsigned long)(nx-4)*ny*sizeof(float));

	// 5âÊëfï™ÇÃïŒç∑ÇÃåvéZ
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx-4 ; j++) {
			max = min = img[j];
			for(k = 1 ; k < 5 ; k++) {
				if(img[i*nx+j+k] > max)      max = img[i*nx+j+k];
				else if(img[i*nx+j+k] < min) min = img[i*nx+j+k];
			}
			dif[i*(nx-4)+j] = max-min;
		}
	}

	// ç≈ëÂílÇ∆ç≈è¨ílÇÃéZèo
	max = min = dif[0];
	for(i = 1 ; i < (nx-4)*ny ; i++) {
		if(dif[i] > max)      max = dif[i];
		else if(dif[i] < min) min = dif[i];
	}
	// ãœàÍê´ÇÃåvéZ
	uni = (max-min)/(max+min)*100;
	printf("\n Differential uniformity = %f\n", uni);
	printf("   [max = %f, min = %f]\n", max, min);
	free(dif);
}
