/*  P5-03backproj.c  */

#include  <math.h>
#define   PI  3.14159265358979

void BackProjection(int pp, float *img, int nx, int ny, double plx, double ply, float *prj, int px, int pa, double pl)
// 逆投影を行う関数
// int    pp;   逆投影をPIで行うか2*PIで行うか（1 or 2）
// float  *img; 再構成した画像データ
// int    nx;   画像のマトリクスサイズ（x方向）
// int    ny;   画像のマトリクスサイズ（y方向）
// double plx;  画像のピクセル実長（x方向：cm）
// double ply;  画像のピクセル実長（y方向：cm）
// float  *prj; 投影データ
// int    px;   投影データの動径方向のデータ数
// int    pa;   投影データの角度方向のデータ数
// double pl    投影データの動径方向のピクセル実長（cm）
{
	int     i, j, k, ix;
	double  x0, cx, cy, th, tx, ty, t1, t2;
	float   *bp2;

	for(i = 0 ; i < nx*ny; i++)
		img[i] = 0;
	for(k = 0 ; k < pa ; k++) {
		th = pp*k*PI/pa;
		cx =  cos(th)*plx/pl;
		cy = -sin(th)*ply/pl;
		x0 = -cx*nx/2-cy*ny/2+px/2;
		bp2 = prj+k*px;
		for(i = 0, ty = x0 ; i < ny ; i++, ty += cy) {
			for(j = 0, tx = ty ; j < nx ; j++, tx += cx) {
				ix = (int)tx;
				if(ix < 0 || ix > px-2)     continue;
				t1 = tx-ix;
				t2 = 1-t1;
				img[i*nx+j] += (float)(t1*bp2[ix+1]+t2*bp2[ix]);
			}
		}
	}
	for(i = 0 ; i < nx*ny ; i++)
		img[i] /= pa;
}
