/* P3-08mkphantom.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  7     /* number of parameters + 1 */

typedef struct phan_data {   /* Phantom data */
	double  x0;     /*  X Coordinate */
	double  y0;     /*  Y Coordinate */
	double  a;      /*  Minor Axis */
	double  b;      /*  Major Axis */
	double  ph;     /*  Rotation angle */
	double  d;      /*  Density */
	struct  phan_data  *next; /*  next self pointer */
} PH_DATA;

typedef struct {
	char     f1[50]; /* input parameter file name */
	char     f2[50]; /* output image file name */
	float    *img;   /* ellipse image matrix */
	int      nx;     /* number of width */
	int      ny;     /* number of height */
	double   pl;     /* pixel length (cm/pixel) */
	double   aw;     /* image area width (cm) */
	PH_DATA  *pd;    /* pointer of Phantom data */
} Param;

char *menu[PN] = {
	"Make Phantom image",
	"Input  parameter file name <.pmt>  ",
	"Output image     file name <float> ",
	"Number of width   ",
	"Number of height  ",
	"Pixel length      ",
	"Image area width  ",
};

void read_phantom_data(char *, PH_DATA *);
void write_data(char *, float *, int);
void mkellipse_phantom(float *, int, int, double, double, double, double, double, double);

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
	sprintf( pm->f2, "shepp.img");
	pm->nx = 128;
	pm->ny = 128;
	pm->pl = 0.15625;
	pm->aw = 20.0;

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
		fprintf( stdout, " %s [%f] :", menu[i++], pm->pl );
		if(*gets(dat) != '\0')  pm->pl = atof(dat);
		fprintf( stdout, " %s [%f] :", menu[i++], pm->aw );
		if(*gets(dat) != '\0')  pm->aw = atof(dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
		if((argc--) > 1) pm->nx = atoi( argv[i++] );
		if((argc--) > 1) pm->ny = atoi( argv[i++] );
		if((argc--) > 1) pm->pl = atof( argv[i++] );
		if((argc--) > 1) pm->aw = atof( argv[i++] );
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

	pm->img = (float *)malloc((unsigned long)pm->nx*pm->ny*sizeof(float));

	printf(" *** Read Phantom data   ***\n");
	pm->pd = (PH_DATA *)malloc(sizeof(PH_DATA));
	pm->pd->next = NULL;
	read_phantom_data(pm->f1, pm->pd);

	printf(" *** Making Phantom Image ***\n");
	for(i = 0 ; i < pm->nx*pm->ny ; i++)
		pm->img[i] = 0;
	for( now = pm->pd ; now != NULL ; now = now->next) {
		now->x0 *= pm->aw/2/pm->pl;  // ’·‚³‚ð‰æ‘f‚ÉŠ·ŽZ
		now->y0 *= pm->aw/2/pm->pl;
		now->a  *= pm->aw/2/pm->pl;
		now->b  *= pm->aw/2/pm->pl;
		mkellipse_phantom(pm->img, pm->nx, pm->ny, now->x0, now->y0, now->a, now->b, now->ph, now->d);
	}

	printf(" *** Write Image data   ***\n");
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
}

void read_phantom_data(char *fi, PH_DATA *now)
{
	int    i, k, flag;
	char   dat[256];
	double w[6];
	FILE   *fp;

	/* open Phantom parameter file */
	if((fp = fopen(fi, "r")) == NULL) {
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
	if((fp = fopen(fi, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(img, sizeof(float), size, fp);
	fclose(fp);
}
