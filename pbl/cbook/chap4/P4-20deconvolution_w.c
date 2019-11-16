/* P4-20deconvolution_w.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  6     /* number of parameters + 1 */
#define  PI  3.14159265358979

typedef struct {
	char   f1[50]; /* input  image file name */
	char   f2[50]; /* input  psf   file name */
	char   f3[50]; /* output image file name */
	float  *img;   /* image real data */
	float  *imi;   /* image imaginary data */
	float  *psf;   /* psf   real data */
	float  *psi;   /* psf   imaginary data */
	int    nx;     /* number of width  */
	int    ny;     /* number of height */
} Param;

char *menu[PN] = {
	"2-D deconvolution",
	"Input  image file name <float> ",
	"Input  psf   file name <float> ",
	"Output image file name <float> ",
	"Number of width          ",
	"Number of height         ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void fft2d(int, float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void inverse(float *, float *, int);
void window(float *, float *, int, int);
void deconv(float *, float *, float *, float *, int);

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
	sprintf( pm->f2, "psf.img");
	sprintf( pm->f3, "n1.img");
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
	int     i;
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->imi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->psf = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->psi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->img, pm->nx*pm->ny);
	read_data(pm->f2, pm->psf, pm->nx*pm->ny);
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->imi[i] = pm->psi[i] = 0;

	printf(" *** %s ***\n", menu[0]);
	fft2d(1, pm->img, pm->imi, pm->nx, pm->ny);
	fft2d(1, pm->psf, pm->psi, pm->nx, pm->ny);

	// psf‚Ì‹t”ŒvZ
	inverse(pm->psf, pm->psi, pm->nx*pm->ny);
	// ƒnƒjƒ“ƒO‘‹ŠÖ”
	window(pm->psf, pm->psi, pm->nx, pm->ny);
	write_data("d1.img", pm->psf, pm->nx*pm->ny);

	// ‰æ‘œ‚Æpsf‚Ì‹t”‚ÌŠ|‚¯Z
	deconv(pm->img, pm->imi, pm->psf, pm->psi, pm->nx*pm->ny);

	fft2d(-1, pm->img, pm->imi, pm->nx, pm->ny);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f3, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm->imi);
	free(pm->psf);
	free(pm->psi);
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

void inverse(float *pre, float *pim, int n)
// ‹t”‚ÌŒvZ
{
	int     i;
	double  pa;

	for(i = 0 ; i < n ; i++) {
		if((pa = pre[i]*pre[i]+pim[i]*pim[i]) > 1e-10) {
			pre[i] = (float)(pre[i]/pa);
			pim[i] = (float)(-pim[i]/pa);
		}
		else {
			pre[i] = pim[i] = 1;
		}
	}
}

void window(float *pre, float *pim, int nx, int ny)
// ƒnƒjƒ“ƒO‘‹ŠÖ”i‚ü”g”‚Ì’[‚Ì’l‚ğ’²ß‚µ‚Ä‚¢‚é: 7/16j
{
	int     i, j;
	double  u, v, w, win;

	for(i = 0 ; i < ny ; i++) {
		v = ny/2-i;
		for(j = 0 ; j < nx ; j++) {
			u = j-nx/2;
			w = sqrt(u*u+v*v);
			if(w > nx/2*7./16.) {
				win = 0;
			}
			else {
				win = 0.5*(1.+cos(PI*w/nx*16./7.));
			}
			pre[i*nx+j] *= (float)win;
			pim[i*nx+j] *= (float)win;
		}
	}
}

void deconv(float *fre, float *fim, float *hre, float *him, int n)
// •¡‘f”‚ÌŠ|‚¯Z
{
	int   i;
	float gre, gim;

	for(i = 0 ; i < n ; i++) {
		gre = fre[i]*hre[i]-fim[i]*him[i];
		gim = fre[i]*him[i]+fim[i]*hre[i];
		fre[i] = gre;
		fim[i] = gim;
	}
}
