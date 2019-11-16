/*  P5-02fourier.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  7
#define  PI  3.14159265358979

typedef struct { // 入力変数
	char   f1[50]; /* input file name */
	float *prj;    /* projection data */
	float *prr;    /* zero padded projection data (real) */
	float *pri;    /* zero padded projection data (imaginary) */
	int    px;     /* number of bins (X) */
	int    pa;     /* number of projections (Thita) */
	int    zx;     /* zero-padding px*4 */
	char   f2[50]; /* output file name */
	float *img;    /* reconstructed image data (real) */
	float *imi;    /* reconstructed image data (imaginary) */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
} Param;

char *menu[PN] = { // 入力の際のコメント（入力変数とリンク）
	"Fourier Method reconstruction program",
	"Input  file name <float> ",
	"  Number of bins         ",
	"  Number of projections  ",
	"Output file name <float> ",
	"  Number of matrix  (x)  ",
	"  Number of matrix  (y)  ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void zero_padding(float *, int, int, float *, int, int);
void Polar_Rect(float *, int, int, float *, int, int);
void FFT_Ramp(float *, float *, int, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void fft2d(int, float *, float *, int, int);

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
	sprintf( pm->f1, "n0.prj");
	pm->px = 128;
	pm->pa = 128;
	pm->zx = 512;
	sprintf( pm->f2, "n0.img");
	pm->nx = 128;
	pm->ny = 128;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->px );
		if(*gets(dat) != '\0')  pm->px = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->pa );
		if(*gets(dat) != '\0')  pm->pa = atoi(dat);
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
		if((argc--) > 1) pm->px = atoi( argv[i++] );
		if((argc--) > 1) pm->pa = atoi( argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
	}
	else {
		usage(argc, argv);
	}
	pm->zx = pm->px*4;
}

main(int argc, char *argv[] )
{
	int     i;
	Param   *pm;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->prj = (float *)malloc((unsigned long)pm->px*pm->pa*sizeof(float));
	pm->prr = (float *)malloc((unsigned long)pm->zx*pm->pa*sizeof(float));
	pm->pri = (float *)malloc((unsigned long)pm->zx*pm->pa*sizeof(float));
	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));
	pm->imi = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Image data   ***\n");
	read_data(pm->f1, pm->prj, pm->px*pm->pa);

	printf(" *** Zero-Padding (x4) ***\n");
	zero_padding(pm->prr, pm->zx, pm->pa, pm->prj, pm->px, pm->pa);
	for(i = 0 ; i < pm->zx*pm->pa ; i++)
		pm->pri[i] = 0;

	printf(" *** FFT (Ramp direction) ***\n");
	FFT_Ramp(pm->prr, pm->pri, pm->zx, pm->pa);

	printf(" *** Polar -> Rectangular ***\n");
	Polar_Rect(pm->img, pm->nx, pm->ny, pm->prr, pm->zx, pm->pa);
	Polar_Rect(pm->imi, pm->nx, pm->ny, pm->pri, pm->zx, pm->pa);

	printf(" *** 2D-IFT ***\n");
	fft2d(-1, pm->img, pm->imi, pm->nx, pm->ny);

	printf(" *** Write Image data  ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->prj);
	free(pm->prr);
	free(pm->pri);
	free(pm->img);
	free(pm->imi);
	free(pm);
}

void read_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and read data */
	if(NULL == (fp = fopen(fi, "rb"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(prj, sizeof(float), size, fp);
	fclose(fp);
}

void write_data(char *fi, float *prj, int size)
{
	FILE   *fp;

	/* open file and write data */
	if(NULL == (fp = fopen(fi, "wb"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(prj, sizeof(float), size, fp);
	fclose(fp);
}

void FFT_Ramp(float *prr, float *pri, int zx, int za)
// 動径方向にフーリエ変換する
// float *prr;   投影データの実部  prr[zx*za]
// float *pri;   投影データの虚部  pri[zx*za]
// int    zx;    投影データの動径方向の数
// int    za;    投影データの角度方向の数
{
	float *xr, *xi, *si, *co;
	unsigned short *br;
	int i, j;

	xr = (float *)malloc((unsigned long)zx*sizeof(float));
	xi = (float *)malloc((unsigned long)zx*sizeof(float));
	si = (float *)malloc((unsigned long)zx*sizeof(float)/2);
	co = (float *)malloc((unsigned long)zx*sizeof(float)/2);
	br = (unsigned short *)malloc((unsigned long)zx*sizeof(unsigned short));
	FFTInit(zx, si, co, br);
	for(i = 0 ; i < za ; i++) {
		for(j = 0 ; j < zx/2 ; j++) {    // FFT用のデータの入れ替え
			xr[j] = prr[i*zx+j+zx/2];
			xr[j+zx/2] = prr[i*zx+j];
			xi[j] = xi[j+zx/2] = 0;
		}
		FFT(1, zx, xr, xi, si, co, br);  // 1次元フーリエ変換
		for(j = 0 ; j < zx/2 ; j++) {    // FFT用のデータの入れ替え
			prr[i*zx+j] = xr[j+zx/2];
			prr[i*zx+j+zx/2] = xr[j];
			pri[i*zx+j] = xi[j+zx/2];
			pri[i*zx+j+zx/2] = xi[j];
		}
	}
	free(xr);
	free(xi);
	free(si);
	free(co);
	free(br);
}

void Polar_Rect(float *img, int nx, int ny, float *prr, int zx, int pa)
// 極座標から直交座標に変換する
// float *img;   直交座標の画像データ  img[nx*ny]
// int   nx;     x方向の数
// int   ny;     y方向の数
// float *prr;   極座標の画像データ    prr[zx*pa]
// int   zx;     動径方向の数
// int   pa;     角度方向の数(360度)
{
	int     i, j, xi, ti, ti1;
	double  x, y, xx, th, dx0, dx1, dt0, dt1;

	for(i = 0 ; i < ny; i++) {
		y = (ny/2-i)*zx/(double)ny;
		for(j = 0 ; j < nx ; j++) {
			x = (j-nx/2)*zx/(double)nx;
			xx = sqrt(x*x+y*y)+zx/2;
			xi = (int)xx;
			if(xi < 0 || xi >= zx-1) {
				img[i*nx+j] = 0;
				continue;
			}
			dx0 = xx-xi;
			dx1 = 1-dx0;
			th = atan2(y, x)*pa/(2*PI);
			if(th < 0.)  th += pa;
			ti = (int)th;
			dt0 = th-ti;
			dt1 = 1-dt0;
			ti1 = ti == pa-1 ? 0 : ti+1;
			img[i*nx+j] = (float)(dt0*dx0*prr[(ti1)*zx+xi+1]
			             +dt0*dx1*prr[(ti1)*zx+xi  ]
			             +dt1*dx0*prr[(ti  )*zx+xi+1]
			             +dt1*dx1*prr[(ti  )*zx+xi  ]);
		}
	}
}

void zero_padding(float *prz, int zx, int za, float *prj, int px, int pa)
// プロジェクションの動径方向にゼロを付加する
// float *prz;  ゼロ付加後の投影データ  prz[zx*za]
// int    zx;   投影データの動径方向の数
// int    za;   投影データの角度方向の数
// float *prj;  ゼロ付加前の投影データ  prj[px*pa]
// int    px;   投影データの動径方向の数
// int    pa;   投影データの角度方向の数
// (zx == 4*px)
// (za == pa)
{
	int i, j;

	for(i = 0 ; i < zx*za ; i++)
		prz[i] = 0;

	for(i = 0 ; i < pa ; i++) {
		for(j = 0 ; j < px ; j++) {
			prz[i*zx+zx/2+j-px/2] = prj[i*px+j];
		}
	}
}
