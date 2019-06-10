/* P3-23interpolation.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  6     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input file name */
	char   f2[50]; /* output file name */
	float  *fb;    /* input data */
	float  *fa;    /* interpolation data */
	int    nx;     /* number of data */
	int    na;     /* number of data after interpolation */
	int    in;     /* type of interpolation */
} Param;

char *menu[PN] = {
	"1-D interpolation",
	"Input  file name <txt> ",
	"Output file name <txt> ",
	"Number of data  (before) ",
	"Number of data  (after)  ",
	"Type of interpolation (0:near, 1:linear, 2:3d-poli) ",
};

void read_data_txt(char *, float *, int);
void write_data_txt(char *, float *, int);
void interpolation(float *, int, float *, int, int);

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
	sprintf( pm->f2, "n1.txt");
	pm->nx = 5;
	pm->na = 41;
	pm->in = 0;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->na );
		if(*gets(dat) != '\0')  pm->na = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->in );
		if(*gets(dat) != '\0')  pm->in = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->na = atoi( argv[i++] );
		if((argc--) > 1) pm->in = atoi( argv[i++] );
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

	pm->fb = (float *)malloc((unsigned long)pm->nx*sizeof(float));
	pm->fa = (float *)malloc((unsigned long)pm->na*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data_txt(pm->f1, pm->fb, pm->nx);

	printf(" *** %s ***\n", menu[0]);
	interpolation(pm->fa, pm->na, pm->fb, pm->nx, pm->in);

	printf(" *** Write Image data   ***\n");
	write_data_txt(pm->f2, pm->fa, pm->na);

	free(pm->fb);
	free(pm->fa);
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

void write_data_txt(char *fi, float *img, int size)
{
	int    i;
	FILE   *fp;

	/* open file and write data */
	if((fp = fopen(fi, "w")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	for(i = 0 ; i < size ; i++) {
		fprintf(fp, "%f\n", img[i]);
	}
	fclose(fp);
}

void interpolation(float *fa, int na, float *fb, int nx, int in)
{
	int    i, ix, xm1, xp0, xp1, xp2;
	double x, dx, b0, b1, b2, b3;

	for(i = 0 ; i < na ; i++)
		fa[i] = 0;

	switch(in) {
		case  0: // Å‹ßÚ•âŠÔ
			for(i = 0 ; i < na ; i++) {
				fa[i] = fb[(int)(i*(nx-1)/(double)(na-1)+0.5)];
			}
			break;
		case  1: // üŒ`•âŠÔ
			for(i = 0 ; i < na-1 ; i++) {
				x = i*(nx-1)/(double)(na-1);
				ix = (int)x;
				dx = x-ix;
				fa[i] = (float)((1-dx)*fb[ix]+dx*fb[ix+1]);
			}
			fa[na-1] = fb[nx-1];
			break;
		case 2: // 3ŽŸ‘½€Ž®•âŠÔ
			for(i = 0 ; i < na-1 ; i++) {
				x = i*(nx-1)/(double)(na-1);
				ix = (int)x;
				dx = x-ix;
				xm1 = ix-1<0?0:ix-1;
				xp0 = ix;
				xp1 = ix+1;
				xp2 = ix+2>nx-1?nx-1:ix+2;
				b0 = -dx*(dx-1)*(dx-2)/6;
				b1 = (dx+1)*(dx-1)*(dx-2)/2;
				b2 = -dx*(dx+1)*(dx-2)/2;
				b3 = dx*(dx+1)*(dx-1)/6;
				fa[i] = (float)(b0*fb[xm1]+b1*fb[xp0]+b2*fb[xp1]+b3*fb[xp2]);
			}
			fa[na-1] = fb[nx-1];
	}
}
