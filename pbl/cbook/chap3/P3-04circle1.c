/* P3-04circle1.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	char    fi[50];
	float   img[128*128];
	int     nx = 128, ny = 128, r = 32;
	int     i, j, x, y;
	FILE    *fp;

	// 画像の初期化
	for (i = 0 ; i < nx*ny ; i++)
		img[i] = 0;

	// 円画像の作成
	for (i = 0 ; i < ny ; i++) {
		y = ny/2-i;  // y座標への変換
		for(j = 0 ; j < nx ; j++) {
			x = j-nx/2;  // x座標への変換
			if(x*x+y*y <= r*r)
				img[i*nx+j] = 100;
		}
	}

	// 画像の書き出し
	printf( "Input new file name: " );
	scanf( "%s", fi );
	if ((fp = fopen ( fi, "wb")) == NULL) {
		printf("Error: file open [%s].\n", fi);
		exit (1);
	}
	fwrite(img, sizeof(float), nx*ny, fp);
	fclose (fp);
}
