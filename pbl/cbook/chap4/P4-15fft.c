/*  P4-15fft.c  */

/* --- プログラムの説明 ---
　1次元フーリエ変換用の関数
このプログラムを使うときは5つの配列を用意しておく
xr[nx]    : 1次元フーリエ変換用の実部データ
xi[nx]    : 1次元フーリエ変換用の虚部データ
si[nx/2]  : フーリエ変換で使うサインデータ
co[nx/2]  : フーリエ変換で使うコサインデータ
brv[nx]   : フーリエ変換で使うバタフライ演算データ

nx        : 1次元フーリエ変換のデータ数（2のベキ乗）

関数：
void FFTInit(int nx, float *si, float *co, unsigned short *brv)
  フーリエ変換に必要なサインとコサインとバタフライ演算データを作成する関数
  FFT関数を使う前に1度使用して必要なデータを作成する

void FFT(int ir, int nx, float *xr, float *xi, float *si, float *co, unsigned short *brv)
  1次元フーリエ変換を実行する関数

void fft2d(int ir, float *fr, float *fi, int nx, int ny)
  2次元フーリエ変換を実行する関数
*/

#include  <stdlib.h>
#include  <math.h>
#define   PI  3.14159265358979

void bitrev(int nx, float *xr, float *xi, unsigned short *brv)
// バタフライ演算の入れ替え
// int   nx;   データ数
// float *xr;  実部のデータ  xr[nx]
// float *xi;  虚部のデータ  xi[nx]
// unsigned short *brv;  交換用のデータ  brv[nx]
{
	int     i, j;
	float   a, b;
	for(i = 0; i < nx; i++){
	    j = brv[i];
	    if(i < j){
	        a = xr[i];
	        b = xi[i];
	        xr[i] = xr[j];
	        xi[i] = xi[j];
	        xr[j] = a;
	        xi[j] = b;
	    }
	}
}


void FFT(int ir, int nx, float *xr, float *xi, float *si, float *co, unsigned short *brv)
// 1次元フーリエ変換
// int    ir;  順変換(1)と逆変換(-1)
// int    nx;  1次元FFTのデータ数
// float *xr;  実部のデータ          xr[nx]
// float *xi;  虚部のデータ          xi[nx]
// float *si;  FFT用のサインデータ   si[nx/2]
// float *co;  FFT用のコサインデータ co[nx/2]
// unsigned short *brv;  FFT用の入れ替えデータ brv[nx]
{
	int     i, j, n1, n2=nx, j3, j4, k, l, ll, d=1, g;
	float   a, b, c, s;

	for(l = 1; l <= nx/2; l *= 2, d += d) {
	    g  = 0;
	    ll = n2;
	    n2 /= 2;
	    for(k = 1; k <= n2; k++) {
	        n1 = k-ll;
	        c  = co[g];
	        s  = -ir*si[g];
	        g += d;
	        for(j = ll; j <= nx; j += ll) {
	            j3 = j+n1-1;
	            j4 = j3+n2;
	            a  = xr[j3]-xr[j4];
	            b  = xi[j3]-xi[j4];
	            xr[j3] += xr[j4];  xi[j3] += xi[j4];
	            xr[j4] =  c*a+s*b; xi[j4] =  c*b-s*a;
	        }
	    }
	}

	bitrev(nx, xr, xi, brv);
	if(ir == -1)
	    for(i = 0; i < nx; i++) {
	        xr[i] /= nx;
	        xi[i] /= nx;
	    }
}

int br(int nx, unsigned nn)
// 交換データ作成用関数
// int      nx;  データ数
// unsigned nn;  交換前のデータ番号
{
	unsigned        r,c;
	r=0;
	for(c = 1; c <= (unsigned)nx/2; c <<= 1) {
	    r <<= 1;
	    if((nn&c) != 0)
	         r++;
	}
	return(r);
}

void FFTInit(int nx, float *si, float *co, unsigned short *brv)
// FFT用のデータ作成用の関数
// int    nx;  FFTのデータ数
// float *si;  サインデータ用配列   si[nx/2]
// float *co;  コサインデータ用配列 co[nx/2]
// unsigned short *brv;  交換データ用配列 brv[nx]
{
	double  d=2.0*PI/nx;
	int     i;
	int     br(int, unsigned);
	for(i = 0; i < nx/4; i++) {
	    si[i] = (float)sin(d*i);
	    co[i+nx/4] = -si[i];
	}
	for(i = nx/4; i < nx/2; i++) {
	    si[i] = (float)sin(d*i);
	    co[i-nx/4] = si[i];
	}
	for(i = 0; i < nx; i++)
	    brv[i] = br(nx, (unsigned)i);
}

void fft2d(int ir, float *fr, float *fi, int nx, int ny)
// 2次元フーリエ変換
// int    ir;   順変換(1)と逆変換(-1)
// float *fr;   2次元FFTの実部のデータ fr[nx*ny]
// float *fi;   2次元FFTの虚部のデータ fi[nx*ny]
// int    nx;   x方向のデータ数
// int    ny;   y方向のデータ数
{
	int   i, j;
	float *gr, *gi, *si, *co;
	unsigned short *br;

	// x方向のフーリエ変換
	gr = (float *)malloc((unsigned long)nx*sizeof(float));
	gi = (float *)malloc((unsigned long)nx*sizeof(float));
	si = (float *)malloc((unsigned long)nx/2*sizeof(float));
	co = (float *)malloc((unsigned long)nx/2*sizeof(float));
	br = (unsigned short *)malloc((unsigned long)nx*sizeof(unsigned short));
	FFTInit(nx, si, co, br);
	for(i = 0 ; i < ny ; i++) {
		for(j = 0 ; j < nx/2 ; j++) {
			gr[j] = fr[i*nx+j+nx/2];
			gr[j+nx/2] = fr[i*nx+j];
			gi[j] = fi[i*nx+j+nx/2];
			gi[j+nx/2] = fi[i*nx+j];
		}
		FFT(ir, nx, gr, gi, si, co, br);
		for(j = 0 ; j < nx/2 ; j++) {
			fr[i*nx+j+nx/2] = gr[j];
			fr[i*nx+j] = gr[j+nx/2];
			fi[i*nx+j+nx/2] = gi[j];
			fi[i*nx+j] = gi[j+nx/2];
		}
	}
	free(gr);
	free(gi);
	free(si);
	free(co);
	free(br);

	// y方向のフーリエ変換
	gr = (float *)malloc((unsigned long)ny*sizeof(float));
	gi = (float *)malloc((unsigned long)ny*sizeof(float));
	si = (float *)malloc((unsigned long)ny/2*sizeof(float));
	co = (float *)malloc((unsigned long)ny/2*sizeof(float));
	br = (unsigned short *)malloc((unsigned long)ny*sizeof(unsigned short));
	FFTInit(ny, si, co, br);
	for(j = 0 ; j < nx ; j++) {
		for(i = 0 ; i < ny/2 ; i++) {
			gr[i] = fr[(i+ny/2)*nx+j];
			gr[i+ny/2] = fr[i*nx+j];
			gi[i] = fi[(i+ny/2)*nx+j];
			gi[i+ny/2] = fi[i*nx+j];
		}
		FFT(ir, ny, gr, gi, si, co, br);
		for(i = 0 ; i < ny/2 ; i++) {
			fr[(i+ny/2)*nx+j] = gr[i];
			fr[i*nx+j] = gr[i+ny/2];
			fi[(i+ny/2)*nx+j] = gi[i];
			fi[i*nx+j] = gi[i+ny/2];
		}
	}
	free(gr);
	free(gi);
	free(si);
	free(co);
	free(br);
}
