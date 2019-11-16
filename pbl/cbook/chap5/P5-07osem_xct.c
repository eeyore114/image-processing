/*  P5-07osem_xct.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  11
#define  NI  3
#define  PI  3.14159265358979

typedef struct { // ���͕ϐ�
	char   f1[50]; /* input file name */
	float  *prj;   /* projection data */
	int    px;     /* number of bins (X) */
	int    pa;     /* number of projections (Thita) */
	double pl;     /* Pixel length of bins */
	char   f2[50]; /* output file name */
	float  *img;   /* reconstructed image data */
	int    nx;     /* number of matrix (x) */
	int    ny;     /* number of matrix (y) */
	double plm;    /* Pixel length of matrix */
	int    nit;    /* number of iteration */
	int    subset; /* subset (OSEM) */
} Param;

char *menu[PN] = { // ���͂̍ۂ̃R�����g�i���͕ϐ��ƃ����N�j
	"OS-ML reconstruction",
	"Projection file name <float>         ",
	"    Number of bins                   ",
	"    Number of projections            ",
	"    Pixel length of projections (cm) ",
	"Image file name <float>              ",
	"    Number of matrix  (x)            ",
	"    Number of matrix  (y)            ",
	"    Pixel length of matrix (cm)      ",
	"Number of iteration                  ",
	"Number of subset                     ",
	};

typedef struct {
	int   x;
	float c[NI];
} CIJ;

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void OSEM(float *, int, int, double, float *, int, int, double, int, int);
void make_cij_xct(CIJ *, int, int, int, int);

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
	sprintf( pm->f2, "n1.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->plm = 0.15625;
	pm->nit = 10;
	pm->subset = 8;

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
		fprintf( stdout, " %s [%f] :", menu[i++], pm->plm );
		if(*gets(dat) != '\0')  pm->plm = atof(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->nit );
		if(*gets(dat) != '\0')  pm->nit = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->subset );
		if(*gets(dat) != '\0')  pm->subset = atoi(dat);
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
		if((argc--) > 1) pm->plm = atof( argv[i++] );
		if((argc--) > 1) pm->nit = atoi( argv[i++] );
		if((argc--) > 1) pm->subset = atoi( argv[i++] );
	}
	else {
		usage(argc, argv);
	}

	// Error
	if(pm->pa%pm->subset != 0) {
		fprintf(stderr, "Error: invalid number of sebset. [angle%%subset==0]\n");
		exit(1);
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

	printf(" *** OSEM reconstruction ***\n");
	OSEM(pm->img, pm->nx, pm->ny, pm->plm, pm->prj, pm->px, pm->pa, pm->pl, pm->nit, pm->subset);

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

void forward_projection(float *aprj, int px, int pa, float *img, int nx, int ny, double lxy, CIJ *c, int sub, int subset)
{
	int    i, j, k;
	float  *p;
	CIJ    *cc;

	for(p = aprj+px*sub, cc = c+nx*ny*sub, i = 0 ; i < pa/subset ; i++, p+=px*subset, cc+=nx*ny*subset) {
		for(j = 0 ; j < px ; j++)
			p[j] = 0;
		for(j = 0 ; j < nx*ny ; j++) {
			for(k = 0 ; k < NI ; k++) {
				p[cc[j].x+k] += (float)(cc[j].c[k]*img[j]*lxy);
			}
		}
	}
}

double osem_1(int nx, int ny, double lxy, float *rprj, int px, int pa, int ij, CIJ *c, int sub, int subset)
{
	int     i, j, k;
	double  cc;
	double  a = 0, r = 0;

	for(i = sub ; i < pa ; i+=subset) {
		j = i*nx*ny+ij;
		for(k = 0 ; k < NI ; k++) {
			cc = c[j].c[k];
			r += cc*rprj[i*px+c[j].x+k];
			a += cc;
		}
	}
	if(a == 0.)  return 0.;
	else         return r/a;
}

void OSEM(float *img, int nx, int ny, double lxy, float *prj, int px, int pa, double lp, int n, int subset)
{
	int     i, j, k, m1, m2, *sub;
	char    fi[50];
	float   *aprj, *rprj, *aimg;
	CIJ     *c;

	aprj = (float *)malloc(px*pa*sizeof(float));
	rprj = (float *)malloc(px*pa*sizeof(float));
	aimg = (float *)malloc(nx*ny*sizeof(float));
	c = (CIJ *)malloc(nx*ny*pa*sizeof(CIJ));

	// �T�u�Z�b�g�isubset�j�̏��Ԃ����肷��
	sub = (int *)malloc(subset*sizeof(int));
	k = 0;
	for(i = 0 ; i < 32 ; i++)
		k += (subset >> i) & 1;
	if(k == 1) {
		m1 = 0;
		sub[m1++] = 0;
		for(i = subset, m2 = 1 ; i > 1 ; i/=2, m2*=2) {
			for(j = 0 ; j <  m2 ; j++)
				sub[m1++] = sub[j]+i/2;
		}
	}
	else {
		for(i = 0 ; i < pa/subset ; i++)
			sub[i] = i;
	}
	printf("\n subset [");
	for(i = 0 ; i < subset ; i++)
		printf(" %d", sub[i]);
	printf(" ]\n");

	// �@ ���o�m��Cij���v�Z����
	printf(" *** Make Cij parameter ***\n");
	make_cij_xct(c, nx, ny, px, pa);

	// ml-em itaration
	// �A �����摜�����肷��
	for(i = 0 ; i < ny*nx ; i++) 
		img[i] = 1;
	sprintf(fi, "osem%03d.img", 0);
	write_data(fi, img, nx*ny);  // �����摜�̏o��

	for(i = 0 ; i < n ; i++) {
		fprintf(stderr, "\r *** OSEM iteration [%2d/%2d]", i+1, n);
		for(k = 0 ; k < subset ; k++) {
			// �B �����摜���瓊�e���v�Z����
			forward_projection(aprj, px, pa, img, nx, ny, lxy, c, sub[k], subset);
			// �C ���e�f�[�^yi�ƁC�B�Ōv�Z�������e�Ƃ̔���v�Z����
			for(j = 0 ; j < px*pa ; j++) {
				if(aprj[j] == 0.)  rprj[j] = 0;
				else               rprj[j] = prj[j]/aprj[j];
			}
			// �D �C�Ōv�Z���ꂽ����t���e����
			// �E �t���e�摜���m���̑��a�ŋK�i������
			for(j = 0 ; j < nx*ny ; j++) {
				if(img[j] == 0.)  aimg[j] = 0;
				else  aimg[j] = (float)osem_1(nx, ny, lxy, rprj, px, pa, j, c, sub[k], subset);
			}
			// �F �t���e�摜�������摜��j(k)�Ɋ|���čX�V�摜��j(k+1)���쐬����
			for(j = 0 ; j < nx*ny ; j++)
				img[j] *= aimg[j];
		}
		sprintf(fi, "osem%02d.prj", i+1);
		write_data(fi, aprj, px*pa);
		sprintf(fi, "osem%02d.prr", i+1);
		write_data(fi, rprj, px*pa);
		sprintf(fi, "osem%02d.rat", i+1);
		write_data(fi, aimg, nx*ny);
		sprintf(fi, "osem%02d.img", i+1);
		write_data(fi, img, nx*ny);
	}
	printf("\n");

	free(aprj);
	free(rprj);
	free(aimg);
	free(c);
}

void make_cij_xct(CIJ *c, int nx, int ny, int px, int pa)
// �P��f�����e�f�[�^�ɓ��e�����l�����߂�֐�
// CIJ   *c;    �P��f�̓��e
// int    nx;   �摜��x�����̐�
// int    ny;   �摜��y�����̐�
// int    px;   ���e�̓��a�����̐�
// int    pa;   ���e�̊p�x�����̐�
{
	int     i, j, k, ix, ij;
	double  x, y, xx, th, a, b, x05, d, si, co;

	for(i = 0 ; i < nx*ny*pa ; i++) {
		c[i].x = 0;
		c[i].c[0] = 0;
		c[i].c[1] = 0;
		c[i].c[2] = 0;
	}

	for(ij = 0, k = 0 ; k < pa ; k++) {
		th = 2*PI*k/pa;
		si = sin(th);
		co = cos(th);
		if(fabs(si) > fabs(co)) {
			a = fabs(si);
			b = fabs(co);
		}
		else {
			a = fabs(co);
			b = fabs(si);
		}
		for(i = 0 ; i < ny ; i++) {
			y = ny/2-i;
			for(j = 0 ; j < nx ; j++, ij++) {
				x = j-nx/2;
				xx = x*co+y*si;
				ix = (int)floor(xx+.5);
				if(ix+nx/2 < 1 || ix+nx/2 > nx-2) continue;
				x05 = ix-.5;
				if((d = x05-(xx-(a-b)/2)) > 0.)
					c[ij].c[0] = (float)(b/(2*a)+d/a);
				else if((d = x05-(xx-(a+b)/2)) > 0.)
					c[ij].c[0] = (float)(d*d/(2*a*b));
				x05 = ix+.5;
				if((d = xx+(a-b)/2-x05) > 0.)
					c[ij].c[2] = (float)(b/(2*a)+d/a);
				else if ((d = xx+(a+b)/2-x05) > 0.)
					c[ij].c[2] = (float)(d*d/(2*a*b));
				c[ij].c[1] = (float)(1.-c[ij].c[0]-c[ij].c[2]);
				c[ij].x = ix+px/2-1;
			}
		}
	}
}
