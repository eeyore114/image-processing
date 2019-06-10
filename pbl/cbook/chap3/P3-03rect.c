/* P3-03rect.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	float   img[128*128];
	int     nx = 128, ny = 128;
	int     i, j;
	FILE    *fp;

	// �摜�̏�����
	for (i = 0 ; i < nx*ny ; i++)
		img[i] = 0;

	// ��`�摜�̍쐬
	for (i = 32 ; i <= 96 ; i++) {
		for(j = 32 ; j <= 96 ; j++) {
			img[i*nx+j] = 1;
		}
	}

	// �摜�̏����o��
	printf( "Input new file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(img, sizeof(float), nx*ny, fp);
	fclose (fp);
}
