/* P3-02ftos.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PN  3     /* number of parameters + 1 */

typedef struct {
	char   f1[50]; /* input image file name */
	char   f2[50]; /* output new image file name */
} Param;

char *menu[PN] = {
	"Convert float image to short image",
	"Input  image     file name <float> ",
	"Output new image file name <short> ",
};

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
	sprintf( pm->f2, "n1.img");

	i = 0;
	if( argc == 1+i ) {
		fprintf( stdout, "\n%s\n\n", menu[i++] );
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f1 );
		if(*gets(dat) != '\0')  strcpy(pm->f1, dat);
		fprintf( stdout, " %s [%s] :", menu[i++], pm->f2 );
		if(*gets(dat) != '\0')  strcpy(pm->f2, dat);
	}
	else if ( argc == PN+i ) {
		fprintf( stderr, "\n%s [%s]\n", argv[i++], menu[0] );
		if((argc--) > 1) strcpy( pm->f1, argv[i++] );
		if((argc--) > 1) strcpy( pm->f2, argv[i++] );
	}
	else {
		usage(argc, argv);
	}

}

main(int argc, char *argv[] )
{
	Param   *pm;
	FILE    *fp1, *fp2;
	float	f;
	short	s;

	pm = (Param *)malloc(sizeof(Param));
	getparameter(argc, argv, pm);

	printf(" *** Read/Write Image data   ***\n");
	/* open file for read data */
	if((fp1 = fopen(pm->f1, "rb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", pm->f1);
		exit(1);
	}
	/* open file for write data */
	if((fp2 = fopen(pm->f2, "wb")) == NULL) {
		fprintf( stderr," Error : file open [%s].\n", pm->f2);
		exit(1);
	}
	// ファイルの最後まで読み出し続ける
	while(fread(&f, sizeof(float), 1, fp1)) { 
		// short型に変換
		f += (float)0.5; // 小数点以下は四捨五入のため
		if(f > 32767.)       s =  32767;
		else if(f < -32768.) s = -32768;
		else                 s =  (short)f;
		                  
		// short型で書き込み
		fwrite(&s, sizeof(short), 1, fp2);
	}
	fclose(fp1);
	fclose(fp2);
	free(pm);
}
