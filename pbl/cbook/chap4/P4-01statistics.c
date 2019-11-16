/* P4-01statistics.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  8

typedef struct {
	char   f1[50]; /* input  file name */
	float  *img;   /* image matrix */
	int    nx;     /* x size (width) */
	int    ny;     /* y size (height) */
	int    x0;     /* x point */
	int    y0;     /* y point */
	int    w;      /* width of calculation */
	int    h;      /* height of calculation */
} Param;

char *menu[PN] = {
	"Make profile for an image",
	"Input  image   file name ",
	"x-Size (width)        ",
	"y-Size (height)       ",
	"x-point               ",
	"y-point               ",
	"width  of calculation ",
	"height of calculation ",
};

void read_data(char *, float *, int);
void statistics(float *, int, int, int, int, int, int);

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
		fprintf( stdout, " %s [%d] :", menu[i++], pm->w );
		if(*gets(dat) != '\0')  pm->w = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->h );
		if(*gets(dat) != '\0')  pm->h = atoi(dat);
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

	printf(" *** Write profile data   ***\n");
	statistics(pm->img, pm->nx, pm->ny, pm->x0, pm->y0, pm->w, pm->h);

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

void statistics(float *img, int nx, int ny, int x0, int y0, int w, int h)
{
	int    i, j, count;
	double total, max, min, average, avedev, sqdev, dev, stdev;

	printf(" (x,y) = (%3d, %3d),", x0, y0);
	printf(" width = %3d, height = %3d\n", w, h);

	total = count = 0;
	max = min = img[y0*nx+x0];
	for(i = y0 ; i < y0+h ; i++) {
		for(j = x0 ; j < x0+w ; j++) {
			total += img[i*nx+j];
			if(max < (double)img[i*nx+j]) max = img[i*nx+j];
			if(min > (double)img[i*nx+j]) min = img[i*nx+j];
			count++;
		}
	}
	average = total/count;
	avedev = sqdev = 0;
	for(i = y0 ; i < y0+h ; i++) {
		for(j = x0 ; j < x0+w ; j++) {
			avedev += fabs(img[i*nx+j]-average);
			sqdev  += (img[i*nx+j]-average)*(img[i*nx+j]-average);
		}
	}
	printf(" number of pixels   = %d\n", count);   // 画素数
	printf(" total counts       = %f\n", total);   // 全カウント数
	printf(" maximum            = %f\n", max);     // 最大値
	printf(" minimum            = %f\n", min);     // 最小値
	printf(" average            = %f\n", average); // 平均値
	printf(" average deviation  = %f\n", avedev);  // 平均偏差
	printf(" square  deviation  = %f\n", sqdev);   // 偏差平方和
	printf(" deviation          = %f\n", dev = sqdev/count); // 分散
	printf(" standard deviation = %f\n", stdev = sqrt(dev)); // 標準偏差（RMSE）
	printf(" %%RMSU              = %f\n", 100*stdev/average); // %RMSU（%RMSE）
}
