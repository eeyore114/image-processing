/*  P5-06mlem_xct.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PN  10
#define  NI  3
#define  PI  3.14159265358979

typedef struct { // “ü—Í•Ï”
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
} Param;

char *menu[PN] = { // “ü—Í‚ÌÛ‚ÌƒRƒƒ“ƒgi“ü—Í•Ï”‚ÆƒŠƒ“ƒNj
	"EM-ML reconstruction",
	"Projection file name <float>         ",
	"    Number of bins                   ",
	"    Number of projections            ",
	"    Pixel length of projections (cm) ",
	"Image file name <float>              ",
	"    Number of matrix  (x)            ",
	"    Number of matrix  (y)            ",
	"    Pixel length of matrix (cm)      ",
	"Number of iteration                  ",
};

typedef struct {
	int   x;
	float c[NI];
} CIJ;

void read_data(char *, float *, int);
void write_data(char *, float *, int);
void ML_EM(float *, int, int, double, float *, int, int, double, int);
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
	pm->nit = 50;

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

	printf(" *** %s ***\n", menu[0]);
	ML_EM(pm->img, pm->nx, pm->ny, pm->plm, pm->prj, pm->px, pm->pa, pm->pl, pm->nit);

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

void forward_projection(float *aprj, int px, int pa, float *img, int nx, int ny, double lxy, CIJ *c)
{
	int    i, j, k;
	float  *p;
	CIJ    *cc;

	for(i = 0 ; i < px*pa ; i++)
		aprj[i] = 0;
	for(p = aprj, cc = c, i = 0 ; i < pa ; i++, p+=px, cc+=nx*ny) {
		for(j = 0 ; j < nx*ny ; j++) {
			for(k = 0 ; k < NI ; k++) {
				p[cc[j].x+k] += (float)(cc[j].c[k]*img[j]*lxy);
			}
		}
	}
}

double ml_em_1(int nx, int ny, double lxy, float *rprj, int px, int pa, int ij, CIJ *c)
{
	int     i, j, k;
	double  cc;
	double  a = 0, r = 0;

	for(i = 0 ; i < pa ; i++) {
		j = i*nx*ny+ij;
		if(c[j].x <= 0 || c[j].x >= px-2)  continue;
		for(k = 0 ; k < NI ; k++) {
			cc = c[j].c[k];
			r += cc*rprj[i*px+c[j].x+k];
			a += cc;
		}
	}
	if(a == 0.)  return 0.;
	else         return r/a;
}

void ML_EM(float *img, int nx, int ny, double lxy, float *prj, int px, int pa, double lp, int n)
{
	int     i, j;
	char    fi[50];
	float   *aprj, *rprj, *aimg;
	CIJ     *c;

	aprj = (float *)malloc(px*pa*sizeof(float));
	rprj = (float *)malloc(px*pa*sizeof(float));
	aimg = (float *)malloc(nx*ny*sizeof(float));
	c = (CIJ *)malloc(nx*ny*pa*sizeof(CIJ));

	// ‡@ ŒŸoŠm—¦Cij‚ğŒvZ‚·‚é
	printf(" *** Make Cij parameter ***\n");
	make_cij_xct(c, nx, ny, px, pa);

	// ml-em itaration
	// ‡A ‰Šú‰æ‘œ‚ğ‰¼’è‚·‚é
	for(i = 0 ; i < ny*nx ; i++) 
		img[i] = 1;
	sprintf(fi, "mlem%03d.img", 0);
	write_data(fi, img, nx*ny);  // ‰Šú‰æ‘œ‚Ìo—Í

	for(i = 0 ; i < n ; i++) {
		fprintf( stderr, "\r *** ML-EM iteration [%2d/%2d]", i+1, n);
		// ‡B ‰Šú‰æ‘œ‚©‚ç“Š‰e‚ğŒvZ‚·‚é
		forward_projection(aprj, px, pa, img, nx, ny, lxy, c);
		// ‡C “Š‰eƒf[ƒ^yi‚ÆC‡B‚ÅŒvZ‚µ‚½“Š‰e‚Æ‚Ì”ä‚ğŒvZ‚·‚é
		for(j = 0 ; j < px*pa ; j++) {
			if(aprj[j] == 0.) rprj[j] = 0;
			else              rprj[j] = prj[j]/aprj[j];
		}
		// ‡D ‡C‚ÅŒvZ‚³‚ê‚½”ä‚ğ‹t“Š‰e‚·‚é
		// ‡E ‹t“Š‰e‰æ‘œ‚ğŠm—¦‚Ì‘˜a‚Å‹KŠi‰»‚·‚é
		for(j = 0 ; j < nx*ny ; j++) {
			if(img[j] == 0.)  aimg[j] = 0;
			else  aimg[j] = (float)ml_em_1(nx, ny, lxy, rprj, px, pa, j, c);
		}
		// ‡F ‹t“Š‰e‰æ‘œ‚ğ‰Šú‰æ‘œƒÉj(k)‚ÉŠ|‚¯‚ÄXV‰æ‘œƒÉj(k+1)‚ğì¬‚·‚é
		for(j = 0 ; j < nx*ny ; j++)
			img[j] *= aimg[j];

		if(i<10 || i%10==9) {  // “r’†‰æ‘œ‚Ìo—Í
			sprintf(fi, "mlem%03d.prj", i+1);
			write_data(fi, aprj, px*pa);
			sprintf(fi, "mlem%03d.prr", i+1);
			write_data(fi, rprj, px*pa);
			sprintf(fi, "mlem%03d.rat", i+1);
			write_data(fi, aimg, nx*ny);
			sprintf(fi, "mlem%03d.img", i+1);
			write_data(fi, img, nx*ny);
		}
	}
	printf("\n");

	free(aprj);
	free(rprj);
	free(aimg);
	free(c);
}

void make_cij_xct(CIJ *c, int nx, int ny, int px, int pa)
// ‚P‰æ‘f‚ª“Š‰eƒf[ƒ^‚É“Š‰e‚³‚ê‚é’l‚ğ‹‚ß‚éŠÖ”
// CIJ   *c;    ‚P‰æ‘f‚Ì“Š‰e
// int    nx;   ‰æ‘œ‚Ìx•ûŒü‚Ì”
// int    ny;   ‰æ‘œ‚Ìy•ûŒü‚Ì”
// int    px;   “Š‰e‚Ì“®Œa•ûŒü‚Ì”
// int    pa;   “Š‰e‚ÌŠp“x•ûŒü‚Ì”
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
