/* P4-22correlation_r.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  6     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input  f(x,y) file name */
	char   f2[50]; /* input  g(x,y) file name */
	char   f3[50]; /* output r(x,y) file name */
	float  *fxy;   /* f(x,y) data */
	float  *gxy;   /* g(x,y) data */
	float  *rxy;   /* r(x,y) data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"2-D correlation (Real space)",
	"Input  f(x,y) file name <float> ",
	"Input  g(x,y) file name <float> ",
	"Output r(x,y) file name <float> ",
	"Number of width   ",
	"Number of height  ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void correlation_r(float *, float *, float *, int, int);

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
	sprintf( pm->f2, "g0.img");
	sprintf( pm->f3, "r0.img");
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
	pm->gxy = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->rxy = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->fxy, pm->nx*pm->ny);
	read_data(pm->f2, pm->gxy, pm->nx*pm->ny);

	printf(" *** %s ***\n", menu[0]);
	correlation_r(pm->rxy, pm->fxy, pm->gxy, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->rxy, pm->nx*pm->ny);

	free(pm->fxy);
	free(pm->gxy);
	free(pm->rxy);
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

void correlation_r(float *rxy, float *fxy, float *gxy, int nx, int ny)
{
	int   i, j, m, n;
	float *gxy2, buff;

	// ëäå›ëää÷ÇÃÇΩÇﬂÅAècâ°2î{Ç…Ç∑ÇÈ
	gxy2 = (float *)malloc((unsigned long)4*nx*ny*sizeof(float));
	// ècâ°2î{Ç…Ç∑ÇÈâÊëúÇ4ï™äÑÇ≈îΩì]Ç≥ÇπÇÈ
	for(i = 0 ; i < ny/2 ; i++) {
		for(j = 0 ; j < nx/2 ; j++) {
			buff = gxy[i*nx+j];
			gxy[i*nx+j] = gxy[(i+ny/2)*nx+j+nx/2];
			gxy[(i+ny/2)*nx+j+nx/2] = buff;
			buff = gxy[i*nx+j+nx/2];
			gxy[i*nx+j+nx/2] = gxy[(i+ny/2)*nx+j];
			gxy[(i+ny/2)*nx+j] = buff;
		}
	}
	// îΩì]âÊëúÇècâ°2î{âÊëúÇÃ4ÉJèäÇ…ìñÇƒÇÕÇﬂÇÈ
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			gxy2[i*nx*2+j] = gxy2[i*nx*2+j+nx] = gxy2[(i+ny)*nx*2+j] = gxy2[(i+ny)*nx*2+j+nx] = gxy[i*nx+j];
		}
	}

	// 2éüå≥ëäå›ëää÷ÇÃé¿çs
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx ; j++) {
			rxy[i*nx+j] = 0;
			for(m = 0 ; m < ny ; m++) {
				for(n = 0 ; n < nx ; n++) {
					rxy[i*nx+j] += fxy[m*nx+n]*gxy2[(i+m)*nx*2+j+n];
				}
			}
		}
	}

	free(gxy2);
}
