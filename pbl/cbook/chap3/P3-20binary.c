/* P3-20binary.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  5     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input image file name */
	char   f2[50]; /* output new image file name */
	float  *img;   /* image matrix */
	int    nx;     /* number of width */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"Ascii to float image",
	"Input  image     file name <ascii> ",
	"Output new image file name <float> ",
	"Number of width  ",
	"Number of height ",
};

void read_data_ascii(char *, float *, int, int);
void write_data(char *, float *, int);

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
	sprintf( pm->f1, "n0.csv");
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
	read_data_ascii(pm->f1, pm->img, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
}

void read_data_ascii(char *fi, float *img, int nx, int ny)
{
	int    i, j, k, r;
	char   c, buff[256];
	FILE   *fp;
	/* open file and write data */
	if((fp = fopen(fi, "r")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	for(i = 0 ; i < nx*ny ; i++)
		img[i] = 0;
	for(i = 0 ; i < ny ; i++) {
		r = 0;
		for(j = 0 ; j < nx ; j++) {
			k = 0;
			while(1) {
				c = getc(fp);
				switch(c) {
					case '\n': case EOF: // x方向の終了判定
						r = 1;
						buff[k] = '\0';
						img[i*nx+j] = (float)atof(buff);
						k = 0;
						break;
					case ',':  // 1画素の終了判定
						buff[k] = '\0';
						img[i*nx+j] = (float)atof(buff);
						k = 0;
						break;
					default:  // 数字データの入力
						buff[k++] = c;
				}
				if(k == 0 || k > 255) break; // 1画素の終了
			}
			if(r == 1) break;  // x方向の終了
		}
	}
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
