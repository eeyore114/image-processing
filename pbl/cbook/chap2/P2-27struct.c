/* P2-27struct.c */

#include  <stdio.h>
#include  <stdlib.h>

typedef struct {
	char    f1[50];
	char    f2[50];
	short	*img;
	int	nx;
	int	ny;
} PARAM;

void read_data(char *, short *, int);
void write_data(char *, short *, int);

void getparameter(PARAM *pm)
{
	char   dat[256];
	printf( "Input image file name : " );
	gets( pm->f1 );
	printf( "Input  new  file name : " );
	gets( pm->f2 );
	printf( "Input image width  : " );
	gets( dat );
	pm->nx = atoi( dat );
	printf( "Input image height : " );
	gets( dat );
	pm->ny = atoi( dat );
}

main( )
{
	PARAM  *pm;

	pm = (PARAM *)malloc(sizeof(PARAM));
	getparameter(pm);

	pm->img = (short *)malloc(pm->nx*pm->ny*sizeof(short));

	read_data(pm->f1, pm->img, pm->nx*pm->ny);
	write_data(pm->f2, pm->img, pm->nx*pm->ny);

	free(pm->img);
	free(pm);
}

void read_data(char *fi, short *buff, int size)
{
	FILE  *fp;

	if ((fp = fopen ( fi, "rb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fread(buff, sizeof(short), size, fp);
	fclose (fp);
}

void write_data(char *fi, short *buff, int size)
{
	FILE  *fp;

	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(buff, sizeof(short), size, fp);
	fclose (fp);
}
