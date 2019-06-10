/* P4-02fwhm_lsf.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  3     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input lsf file name */
	float  *img;   /* lsf data */
	int    nx;     /* number of data */
} Param;

char *menu[PN] = {
	"FWHM and FWTM of lsf",
	"Input lsf file name <txt> ",
	"Number of data            ",
};

void read_data_txt(char *, float *, int);
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
	sprintf( pm->f1, "n0.txt");
	pm->nx = 128;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
	}
	else {
		usage(argc, argv);
	}

}

main(int argc, char *argv[] )
{
	Param   *pm;
	double  fwhm, fwtm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->img = (float *)malloc((unsigned long)pm->nx*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data_txt(pm->f1, pm->img, pm->nx);

	printf(" *** FWHM and FWTM of lsf ***\n");
	fwhm = fwxm_lsf(pm->img, pm->nx, 0.5);
	fwtm = fwxm_lsf(pm->img, pm->nx, 0.1);

	printf(" FWHM = %f\n", fwhm);
	printf(" FWTM = %f\n", fwtm);

	free(pm->img);
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

	for(i = mx ; (double)p[i] > hv && i < nx-1 ; i++);
	dx = (p[i-1] == p[i]? 0.5 : (p[i-1]-hv)/(p[i-1]-p[i]));
	fw1 = i-1+dx; // 最大値から右側の幅

//	fprintf( stderr,"  max*%f=%8.1f,", rt, hv);
//	fprintf( stderr," p[%3d]=%8.1f, p[%3d]=%8.1f, fw1=%8.3f\n", i, p[i-1], i+1, p[i], fw1);

	return fw1-fw0;
}
