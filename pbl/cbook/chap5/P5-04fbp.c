/*  P5-04fbp.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  10
#define  PI  3.14159265358979

typedef struct { // ���͕ϐ�
	char   f1[50]; /* input file name */
	float  *prj;   /* projection data */
	int    px;     /* number of bins (X) */
	int    pa;     /* number of projections (Theta) */
	double pl;     /* Pixel length of bins */
	char   f2[50]; /* output file name */
	float  *img;   /* reconstructed image data */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
	double plx;    /* Pixel length of x-axis */
	double ply;    /* Pixel length of y-axis */ 
} Param;

char *menu[PN] = { // ���͂̍ۂ̃R�����g�i���͕ϐ��ƃ����N�j
	"Filtered Back-Projection",
	"Projection file name <float>       ",
	"  Number of bins                   ",
	"  Number of projections            ",
	"  Pixel length of projections (cm) ",
	"Image file name <float>            ",
	"  Number of matrix  (x)            ",
	"  Number of matrix  (y)            ",
	"  Pixel length of x-axis (cm)      ",
	"  Pixel length of y-axis (cm)      ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void FBP(float *, int, int, double, double, float *, int, int, double);
void zero_padding(float *, int, int, float *, int, int);
void FFT_Filter_IFT(float *, int, int, double);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void Filter(float *, int, double);
void fswap_1D(float *, int);
void BackProjection(int, float *, int, int, double, double, float *, int, int, double);

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
	pm->pl = 0.15625;
	sprintf( pm->f2, "n0.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->plx = 0.15625;
	pm->ply = 0.15625;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->px );
		if(*gets(dat) != '\0')  pm->px = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->pa );
		if(*gets(dat) != '\0')  pm->pa = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->pl );
		if(*gets(dat) != '\0')  pm->pl = atof(dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nx );
		if(*gets(dat) != '\0')  pm->nx = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->ny );
		if(*gets(dat) != '\0')  pm->ny = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->plx );
		if(*gets(dat) != '\0')  pm->plx = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->ply );
		if(*gets(dat) != '\0')  pm->ply = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) pm->px = atoi( argv[i++] );
		if((argc--) > 1) pm->pa = atoi( argv[i++] );
		if((argc--) > 1) pm->pl = atof( argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->plx = atof( argv[i++] );
		if((argc--) > 1) pm->ply = atof( argv[i++] );
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

	pm->prj = (float *)malloc((unsigned long)pm->px*pm->pa*sizeof(float));
	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Projection data   ***\n");
	read_data(pm->f1, pm->prj, pm->px*pm->pa);

	printf(" *** Filtered Back-Projection ***\n");
	FBP(pm->img, pm->nx, pm->ny, pm->plx, pm->ply, pm->prj, pm->px, pm->pa, pm->pl);

	printf(" *** Write Image data  ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->prj);
	free(pm->img);
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

void Filter(float *xr, int nx, double pl)
// Ramp�t�B���^���|����֐�
// float *xr;  �t�B���^���|����f�[�^�z�� xr[nx]
// int    nx;  �t�B���^���|����f�[�^��
// double pl;  �f�[�^��1�s�N�Z���̒���(cm)
{
	int     i;
	double  h;

	h = PI/nx/pl;
	for(i = 0 ; i < nx/2 ; i++)
		xr[i] *= (float)(i*h);
	for(i = nx/2 ; i < nx ; i++)
		xr[i] *= (float)((nx-i)*h);
}

void FFT_Filter_IFT(float *prz, int zx, int za, double pl)
// �t�[���G�̈�œ��e�f�[�^�Ƀt�B���^���|����
// float *prz;   ���e�f�[�^ prz[zx*za]
// int    zx;    ���a�����̐�
// int    za;    �p�x�����̐�
// double pl;    ���a�����̂P�s�N�Z���̒���(cm)
{
	float *xr, *xi, *si, *co;
	unsigned short *br;
	int   i, j;
	void  fswap_1D(float *, int);

	xr = (float *)malloc((unsigned long)zx*sizeof(float));
	xi = (float *)malloc((unsigned long)zx*sizeof(float));
	si = (float *)malloc((unsigned long)zx*sizeof(float)/2);
	co = (float *)malloc((unsigned long)zx*sizeof(float)/2);
	br = (unsigned short *)malloc((unsigned long)zx*sizeof(unsigned short));
	FFTInit(zx, si, co, br);
	for(i = 0 ; i < za ; i++) {
		for(j = 0 ; j < zx/2 ; j++) {    // FFT�p�̃f�[�^�̓���ւ�
			xr[j] = prz[i*zx+j+zx/2];
			xr[j+zx/2] = prz[i*zx+j];
			xi[j] = xi[j+zx/2] = 0;
		}
		FFT(1, zx, xr, xi, si, co, br);  // 1�����t�[���G�ϊ�
		Filter(xr, zx, pl);              // �����f�[�^�ւ̃t�B���^�����O
		Filter(xi, zx, pl);              // �����f�[�^�ւ̃t�B���^�����O
		FFT(-1, zx, xr, xi, si, co, br); // 1�����t�t�[���G�ϊ�
		for(j = 0 ; j < zx/2 ; j++) {    // FFT�p�̃f�[�^�̓���ւ�
			prz[i*zx+j] = xr[j+zx/2];
			prz[i*zx+j+zx/2] = xr[j];
		}
	}
	free(xr);
	free(xi);
	free(si);
	free(co);
	free(br);
}

void FBP(float *img, int nx, int ny, double plx, double ply, float *prj, int px, int pa, double pl)
// �t�B���^�␳�t���e(FBP)�p�̊֐�
// float *img;  �č\�������摜�f�[�^
// int nx;      �摜�̃}�g���N�X�T�C�Y�ix�����j
// int ny;      �摜�̃}�g���N�X�T�C�Y�iy�����j
// double plx;  �摜�̃s�N�Z�������ix�����Fcm�j
// double ply;  �摜�̃s�N�Z�������iy�����Fcm�j
// float *prj;  ���e�f�[�^
// int px;      ���e�f�[�^�̓��a�����̃f�[�^��
// int pa;      ���e�f�[�^�̊p�x�����̃f�[�^��
// double pl;   ���e�f�[�^�̓��a�����̃s�N�Z�������icm�j
{
	int   pz = px*2; // 4�{�[���p�f�B���O�����T���v�����O��
	float *prz;

	prz = (float *)malloc((unsigned long)pz*pa*sizeof(float));

	printf(" *** Zero-Padding (x4) ***\n");
	zero_padding(prz, pz, pa, prj, px, pa);

	printf(" *** FFT -> Filter -> IFT ***\n");
	FFT_Filter_IFT(prz, pz, pa, pl);
	
	printf(" *** Back-Projection   ***\n");
	BackProjection(2, img, nx, ny, plx, ply, prz, pz, pa, pl);

	free(prz);
}

void zero_padding(float *prz, int zx, int za, float *prj, int px, int pa)
// �v���W�F�N�V�����̓��a�����Ƀ[����t������
// float *prz;  �[���t����̓��e�f�[�^  prz[zx*za]
// int    zx;   ���e�f�[�^�̓��a�����̐�
// int    za;   ���e�f�[�^�̊p�x�����̐�
// float *prj;  �[���t���O�̓��e�f�[�^  prj[px*pa]
// int    px;   ���e�f�[�^�̓��a�����̐�
// int    pa;   ���e�f�[�^�̊p�x�����̐�
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
