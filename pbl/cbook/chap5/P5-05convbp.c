/*  P5-05convbp.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  11
#define  PI  3.14159265358979

typedef struct { // ���͕ϐ�
	char   f1[50]; /* input file name */
	float  *prj;   /* projection data */
	int    px;     /* number of bins (X) */
	int    pa;     /* number of projections (Thita) */
	double pl;     /* Pixel length of bins */
	int    nc;     /* number of convolution length */
	char   f2[50]; /* output file name */
	float  *img;   /* reconstructed image data */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
	double plx;    /* Pixel length of x-axis */
	double ply;    /* Pixel length of y-axis */
} Param;

char *menu[PN] = { // ���͂̍ۂ̃R�����g�i���͕ϐ��ƃ����N�j
	"Convolution Back-Projection",
	"Projection file name <float>       ",
	"  Number of bins                   ",
	"  Number of projections            ",
	"  Pixel length of projections (cm) ",
	"Number of convolution length       ",
	"Image file name <float>            ",
	"  Number of matrix  (x)            ",
	"  Number of matrix  (y)            ",
	"  Pixel length of x-axis (cm)      ",
	"  Pixel length of y-axis (cm)      ",
};

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void write_profile(char *, float *, int);
void zero_correct(float *, int, int);
void CBP(float *, int, int, double, double, float *, int, int, double, int);
void period_padding(float *, int, int, float *, int, int);
void Convolution(float *, int, int, double, int);
void FFTInit(int, float *, float *, unsigned short *);
void FFT(int, int, float *, float *, float *, float *, unsigned short *);
void Make_Filter(float *, int, double);
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
	pm->nc = 128;
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
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nc );
		if(*gets(dat) != '\0')  pm->nc = atoi(dat);
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
		if((argc--) > 1) pm->nc = atoi( argv[i++] );
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

	printf(" *** Convolution Back-Projection ***\n");
	CBP(pm->img, pm->nx, pm->ny, pm->plx, pm->ply, pm->prj, pm->px, pm->pa, pm->pl, pm->nc);

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

void CBP(float *img, int nx, int ny, double plx, double ply, float *prj, int px, int pa, double pl, int nc)
// �d���ϕ��ɂ��t�B���^�␳�t���e�@(Convolution Back-Projection)
// float *img;  �č\�������摜�f�[�^
// int    nx;   �摜�̃}�g���b�N�X�T�C�Y�ix�����j
// int    ny;   �摜�̃}�g���b�N�X�T�C�Y�iy�����j
// double plx;  �摜�̃s�N�Z�������ix�����Fcm�j
// double ply;  �摜�̃s�N�Z�������iy�����Fcm�j
// float *prj;  �v���W�F�N�V�����̃f�[�^
// int    px;   �v���W�F�N�V�����̓��a�����̃f�[�^��
// int    pa;   �v���W�F�N�V�����̊p�x�����̃f�[�^��
// double pl;   �v���W�F�N�V�������a�����̃s�N�Z�������icm�j
// int    nc;   �d���ϕ��֐���x�����̃f�[�^��
{
	int   px2 = px*2; // �d���ϕ��p��2�{�����T���v�����O��
	float *pr2;

	pr2 = (float *)malloc((unsigned long)px2*pa*sizeof(float));

	printf(" *** Period-Padding (x2) ***\n");
	period_padding(pr2, px2, pa, prj, px, pa);

	printf(" *** Convolution ***\n");
	Convolution(pr2, px2, pa, pl, nc);

	printf(" *** Back-Projection   ***\n");
	BackProjection(2, img, nx, ny, plx, ply, pr2, px2, pa, pl);

	free(pr2);
}

void period_padding(float *prz, int zx, int za, float *prj, int px, int pa)
// �v���W�F�N�V�����̓��a������2�{�ɂ��ăf�[�^�������I�ɕt������
// float *prz;  �[���t����̓��e�f�[�^  prz[zx*za]
// int    zx;   ���e�f�[�^�̓��a�����̐�
// int    za;   ���e�f�[�^�̊p�x�����̐�
// float *prj;  �[���t���O�̓��e�f�[�^  prj[px*pa]
// int    px;   ���e�f�[�^�̓��a�����̐�
// int    pa;   ���e�f�[�^�̊p�x�����̐�
// (zx == 2*px)
// (za == pa)
{
	int i, j;

	for(i = 0 ; i < zx*za ; i++)
		prz[i] = 0;

	for(i = 0 ; i < pa ; i++) {
		for(j = 0 ; j < px/2 ; j++) {
			prz[i*zx+j] = prz[i*zx+px+j] = prj[i*px+px/2+j];
			prz[i*zx+px/2+j] = prz[i*zx+3*px/2+j] = prj[i*px+j];
		}
	}
}

void Convolution(float *prz, int zx, int za, double pl, int nc)
// ���̈�ŏd���ϕ��ɂ���ăv���W�F�N�V�����Ƀt�B���^���|����
// float *prz;  ���e�f�[�^ prz[zx*za]
// int    zx;   ���e�f�[�^�̓��a�����̐�
// int    za;   ���e�f�[�^�̊p�x�����̐�
// double pl;   �P�s�N�Z���̒��� (cm)
// int    nc;   �R���{�����[�V�����֐��̒���
{
	float  *xr, *xi, *si, *co, buff;
	unsigned short *br;
	int    i, j, k;

	if(nc >= zx/2) {
		fprintf(stderr, "Warning : convolution function is too large.\n");
		nc = zx/2-1;
	}
	xr = (float *)malloc((unsigned long)zx/2*sizeof(float));
	xi = (float *)malloc((unsigned long)zx/2*sizeof(float));
	si = (float *)malloc((unsigned long)zx/2*sizeof(float)/2);
	co = (float *)malloc((unsigned long)zx/2*sizeof(float)/2);
	br = (unsigned short *)malloc((unsigned long)zx/2*sizeof(unsigned short));
	FFTInit(zx/2, si, co, br);
	Make_Filter(xr, zx/2, pl);       // �R���{�����[�V�����֐��̍쐬
	Make_Filter(xi, zx/2, pl);       // �R���{�����[�V�����֐��̍쐬
	FFT(-1, zx/2, xr, xi, si, co, br);  // 1�����t�t�[���G�ϊ�
	for(i = 0 ; i < zx/4 ; i++) {    // FFT�p�̃f�[�^�̓���ւ�
		buff = xr[i];
		xr[i] = xr[i+zx/4];
		xr[i+zx/4] = buff;
	}
	free(si);
	free(co);
	free(br);

	for(i = 0 ; i < za ; i++) {
		for(j = 0 ; j < zx/2 ; j++) {
			xi[j] = 0;
			for(k = 0 ; k < nc ; k++) {
				xi[j] += prz[i*zx+j+zx/4+k-nc/2]*xr[nc/2-k+zx/4];
			}
		}
		for(j = 0 ; j < zx ; j++)
			prz[i*zx+j] = 0;
		for(j = 0 ; j < zx/2 ; j++)
			prz[i*zx+j+zx/4] = xi[j];
	}
	free(xr);
	free(xi);
}

void Make_Filter(float *xr, int nx, double pl)
// Ramachandran�̃t�B���^���쐬����
// float *xr;  �t�B���^��1�����f�[�^
// int    nx;  �f�[�^��
// double pl;  �P�s�N�Z���̒���(cm)
{
	int     i;
	double  h;

	h = PI/nx/pl;
	for(i = 0 ; i < nx/2 ; i++)
		xr[i] = (float)(i*h);
	for(i = nx/2 ; i < nx ; i++)
		xr[i] = (float)((nx-i)*h);
}
