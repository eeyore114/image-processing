/* P4-03fwhm_psf.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  8     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input lsf file name */
	float  *img;   /* lsf data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
	int    x0;     /* start point(x) of area */
	int    y0;     /* start point(y) of area */
	int    w;      /* width of area */
	int    h;      /* height of area */
} Param;

char *menu[PN] = {
	"FWHM and FWTM of psf",
	"Input psf file name <float> ",
	"Number of width        ",
	"Number of height       ",
	"Start point(x) of area ",
	"Start point(x) of area ",
	"width  of area         ",
	"height of area         ",
};

void read_data(char *, float *, int);
void fwxm_psf(float *, int, int, int, int, int, int);
double fwxm_lsf(float *, int, double);

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
	pm->x0 = 0;
	pm->y0 = 0;
	pm->w  = 128;
	pm->h  = 128;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->x0 );
		if(*gets(dat) != '\0')  pm->x0 = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->y0 );
		if(*gets(dat) != '\0')  pm->y0 = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->w  );
		if(*gets(dat) != '\0')  pm->w  = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->h  );
		if(*gets(dat) != '\0')  pm->h  = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->x0 = atoi( argv[i++] );
		if((argc--) > 1) pm->y0 = atoi( argv[i++] );
		if((argc--) > 1) pm->w  = atoi( argv[i++] );
		if((argc--) > 1) pm->h  = atoi( argv[i++] );
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

	printf(" *** FWHM and FWTM of psf ***\n");
	fwxm_psf(pm->img, pm->nx, pm->ny, pm->x0, pm->y0, pm->w, pm->h);

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

void fwxm_psf(float *img, int nx, int ny, int x0, int y0, int w, int h)
{
	int     i, j;
	float   *p;
	double  fwhm, fwtm;

	printf("*** x-direction ***\n");
	p = (float *)malloc((unsigned long)w*sizeof(float));
	for(i = 0 ; i < w ; i++) {
		p[i] = 0;
		for(j = 0 ; j < h ; j++) {
			p[i] += img[(j+y0)*nx+i+x0];
		}
	}
	fwhm = fwxm_lsf(p, w, 0.5);
	fwtm = fwxm_lsf(p, w, 0.1);
	printf("  FWHM = %f\n", fwhm);
	printf("  FWTM = %f\n", fwtm);
	printf("\n");
	free(p);

	printf("*** y-direction ***\n");
	p = (float *)malloc((unsigned long)h*sizeof(float));
	for(i = 0 ; i < h ; i++) {
		p[i] = 0;
		for(j = 0 ; j < w ; j++) {
			p[i] += img[(i+y0)*nx+j+x0];
		}
	}
	fwhm = fwxm_lsf(p, h, 0.5);
	fwtm = fwxm_lsf(p, h, 0.1);
	printf("  FWHM = %f\n", fwhm);
	printf("  FWTM = %f\n", fwtm);
	printf("\n");
	free(p);
}

double fwxm_lsf(float *p, int nx, double rt)
{
	int     i, mx;
	double  max, dx, fw0, fw1, hv;

	max = p[0];
	for(i = 1 ; i < nx ; i++) { // 最大値とその座標の算出
		if((double)p[i] > max) {
			max = p[i];
			mx = i;
		}
	}

	hv = max*rt; // 割合をかけた値の算出

	for(i = mx ; (double)p[i] > hv && i > 0 ; i--);
	dx = (p[i+1] == p[i]? 0.5 : (hv-p[i])/(p[i+1]-p[i]));
	fw0 = i+dx; // 最大値から左側の幅

//	fprintf( stderr,"  max*%f=%8.1f,", rt, hv);
//	fprintf( stderr," p[%3d]=%8.1f, p[%3d]=%8.1f, fw0=%8.3f\n", i, p[i], i+1, p[i+1], fw0);

	for(i = mx ; p[i] > (double)hv && i < nx-1 ; i++);
	dx = (p[i-1] == p[i]? 0.5 : (p[i-1]-hv)/(p[i-1]-p[i]));
	fw1 = i-1+dx; // 最大値から右側の幅

//	fprintf( stderr,"  max*%f=%8.1f,", rt, hv);
//	fprintf( stderr," p[%3d]=%8.1f, p[%3d]=%8.1f, fw1=%8.3f\n", i, p[i-1], i+1, p[i], fw1);

	return fw1-fw0;
}
