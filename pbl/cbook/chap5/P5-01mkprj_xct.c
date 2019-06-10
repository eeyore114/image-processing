/*  P5-01mkprj_xct.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  PI  3.14159265358979
#define  PN  6     /* number of parameters + 1 */

typedef struct phan_data {   /* Phantom data */
	double  x0;     /*  X Coordinate */
	double  y0;     /*  Y Coordinate */
	double  a;      /*  Minor Axis */
	double  b;      /*  Major Axis */
	double  ph;     /*  Rotation angle */
	double  d;      /*  Density */
	struct phan_data  *next; /*  next self pointer */
} PH_DATA;

typedef struct {
	char     f1[50]; /* input parameter file name */
	char     f2[50]; /* output projection file name */
	float    *prj;   /* projction data */
	int      px;     /* number of bins */
	int      pa;     /* number of projections */
	double   pl;     /* pixel length */
	PH_DATA  *pd;  /* pointer of Phantom data */
} Param;

char *menu[PN] = {
	"Make Projection data for X-CT",
	"Input  parameter  file name <.pmt>  ",
	"Output projection file name <float> ",
	"Number of bins           ",
	"Number of projections    ",
	"Pixel length (cm)        ",
};

void read_phantom_data(char *, PH_DATA *);
void write_data(char *, float *, int);
void make_ellipse_projection(float *, int, int, double, double, double, double, double, double, double);

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
	sprintf( pm->f1, "P3-09shepp.pmt");
	sprintf( pm->f2, "n0.prj");
	pm->px = 128;
	pm->pa = 128;
	pm->pl = 0.15625;

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->px );
		if(*gets(dat) != '\0')  pm->px = atoi(dat);
		fprintf( stdout, " %s [%d] :", menu[i++], pm->pa );
		if(*gets(dat) != '\0')  pm->pa = atoi(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->pl );
		if(*gets(dat) != '\0')  pm->pl = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->px = atoi( argv[i++] );
		if((argc--) > 1) pm->pa = atoi( argv[i++] );
		if((argc--) > 1) pm->pl = atof( argv[i++] );
	}
	else {
		usage(argc, argv);
	}

}

main(int argc, char *argv[] )
{
	int     i;
	Param   *pm;
	PH_DATA *now;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	pm->prj = (float *)malloc((unsigned long)pm->px*pm->pa*sizeof(float));

	printf(" *** Read Phantom data   ***\n");
	pm->pd = (PH_DATA *)malloc(sizeof(PH_DATA));
	pm->pd->next = NULL;
	read_phantom_data(pm->f1, pm->pd);

	printf(" *** Making Projection ***\n");
	for(i = 0 ; i < pm->px*pm->pa ; i++)
		pm->prj[i] = 0;
	for(now = pm->pd ; now != NULL ; now = now->next)
		make_ellipse_projection(pm->prj, pm->px, pm->pa, pm->pl, now->x0, now->y0, now->a, now->b, now->ph, now->d);

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->prj, pm->px*pm->pa);

	free(pm->prj);
	free(pm);
}

void read_phantom_data(char *fi, PH_DATA *now)
{
	int     i, k, flag;
	char    dat[256];
	double  w[6];
	FILE    *fp;

	/* open Phantom parameter file */
	if(NULL == (fp = fopen(fi, "r"))) {
		fprintf( stderr, "Error: file open [%s].\n", fi);
		exit(1);
	}

	/* Input Phatom parameters */
	flag = 0;
	while(fgets(dat,100,fp) != NULL) {
		if(*dat=='#'){
			printf("      ");
			printf(dat);
			continue;
		}
		for(i=0;i<6;i++) w[i]=0;
		k = 0;
		for(i=0;i<6;i++){
		  while((dat[k]==' ')||(dat[k]=='\t')) k++;
		  w[i]=atof(dat+k);
		  while((dat[k]!=' ')&&(dat[k]!='\t')) k++;
		}
		if(flag) {
			now->next = (PH_DATA *)malloc(sizeof(PH_DATA));
			now = now->next;
			now->next = NULL;
		}
		now->x0 = w[0];
		now->y0 = w[1];
		now->a  = w[2];
		now->b  = w[3];
		now->ph = w[4];
		now->d  = w[5];
		flag++;
		printf("* %2d *", flag);
		printf("%8.4f,",  now->x0);
		printf("%8.4f,",  now->y0);
		printf("%8.4f,",  now->a);
		printf("%8.4f,",  now->b);
		printf("%8.4f,",  now->ph);
		printf("%8.4f\n", now->d);
	}
	printf("\n");
	fclose(fp);
}

void write_data(char *fi, float *img, int size)
{
	FILE   *fp;

	/* open file and write data */
	if(NULL == (fp = fopen(fi, "wb"))) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(img, sizeof(float), size, fp);
	fclose(fp);
}

void make_ellipse_projection(float *prj, int px, int pa, double pl, double x0, double y0, double a, double b, double phi, double d)
// float  *prj; 作成されるプロジェクションデータ
// int    px;   プロジェクションの動径方向の数
// int    pa;   プロジェクションの角度方向の数（360度）
// double pl;   プロジェクションの動径方向のピクセル実長（cm/pixel）
// double x0;   楕円の中心のx座標（領域の両端を±1.0に規格化してある）
// double y0;   楕円の中心のy座標（同上）
// double a;    楕円のx軸方向の径（同上）
// double b;    楕円のy軸方向の径（同上）
// double phi;  楕円の傾き（度）
// double d;    楕円の濃度
{
	int     i, j;
	double  x1, y1, ph, a2, b2, theta, tp, co, si, x, alpha, beta, ganma, sq;

	x0 *= px/2; // ピクセル値に変換（±1.0 ⇒ ±px/2）
	y0 *= px/2; // ピクセル値に変換
	a  *= px/2; // ピクセル値に変換
	b  *= px/2; // ピクセル値に変換
	ph = PI*phi/180.; // 度からラジアンへ変換
	x1 =  x0*cos(ph)+y0*sin(ph);
	y1 = -x0*sin(ph)+y0*cos(ph);
	a2 = a*a;
	b2 = b*b;
	for(i = 0 ; i < pa ; i++) {
		theta = 2*PI*i/pa;
		tp = theta-ph;
		co = cos(tp);
		si = sin(tp);
		alpha = a2*co*co+b2*si*si;
		for(j = 0 ; j < px ; j++) {
			x = j-px/2;
			beta  = (a2-b2)*co*si*x+b2*si*x1-a2*co*y1;
			ganma = b2*(x*co-x1)*(x*co-x1)+a2*(x*si-y1)*(x*si-y1)-a2*b2;
			sq = beta*beta-alpha*ganma;
		 	if(sq > 0.0){
				prj[i*px+j] += (float)(2*d*pl*sqrt(sq)/alpha);
			}
		}
	}
}
