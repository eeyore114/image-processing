/* P4-12fourier1d.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  7     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input real part file name */
	char   f2[50]; /* input imaginary part file name */
	char   f3[50]; /* output real part file name */
	char   f4[50]; /* output imaginary part file name */
	float  *fr;    /* real part data */
	float  *fi;    /* imaginary part data */
	int    nx;     /* number of data */
	int    ir;     /* forward or inverse */
} Param;

char *menu[PN] = {
	"1-D fourier transform",
	"Input  real      part file name <txt> ",
	"Input  imaginary part file name <txt> ",
	"Output real      part file name <txt> ",
	"Output imaginary part file name <txt> ",
	"Number of data           ",
	"Forward(1) or Inverse(-1)",
};

void read_data_txt(char *, float *, int);
void write_data_txt(char *, float *, int);
void fourier1d(int, float *, float *, int);

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
	sprintf( pm->f1, "n0r.txt");
	sprintf( pm->f2, "n0i.txt");
	sprintf( pm->f3, "n1r.txt");
	sprintf( pm->f4, "n1i.txt");
	pm->nx = 32;
	pm->ir = 1;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f3 );
		if(*gets(dat) != '\0')  strcpy(pm->f3, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f4 );
		if(*gets(dat) != '\0')  strcpy(pm->f4, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ir );
		if(*gets(dat) != '\0')  pm->ir = atoi(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) strcpy( pm->f3, argv[i++] );
		if((argc--) > 1) strcpy( pm->f4, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ir = atoi( argv[i++] );
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

	pm->fr = (float *)malloc((unsigned long)pm->nx*sizeof(float));
	pm->fi = (float *)malloc((unsigned long)pm->nx*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data_txt(pm->f1, pm->fr, pm->nx);
	read_data_txt(pm->f2, pm->fi, pm->nx);

	printf(" *** 1-D fourier transform   ***\n");
	fourier1d(pm->ir, pm->fr, pm->fi, pm->nx);

	printf(" *** Write Image data   ***\n");
	write_data_txt(pm->f3, pm->fr, pm->nx);
	write_data_txt(pm->f4, pm->fi, pm->nx);

	free(pm->fr);
	free(pm->fi);
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

void fourier1d(int ir, float *fr, float *fi, int nx)
{
	int     i, j, n = 1;
	float   *gr, *gi;
	double  u, x;

	gr = (float *)malloc((unsigned long)nx*sizeof(float));
	gi = (float *)malloc((unsigned long)nx*sizeof(float));

	for(i = 0 ; i < nx ; i++) {
		u = i-nx/2;
		gr[i] = gi[i] = 0;
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;
			gr[i] += (float)( fr[j]*cos(2*PI*u*x/nx)+ir*fi[j]*sin(2*PI*u*x/nx));
			gi[i] += (float)(-ir*fr[j]*sin(2*PI*u*x/nx)+fi[j]*cos(2*PI*u*x/nx));
		}
	}

	if(ir == -1)  n = nx; // ‹t•ÏŠ·‚Íƒf[ƒ^”‚ÅŠ„‚é
	for(i = 0 ; i < nx ; i++) {
		fr[i] = gr[i]/n;
		fi[i] = gi[i]/n;
	}

	free(gr);
	free(gi);
}
