/*  P4-15fft.c  */

/* --- �v���O�����̐��� ---
�@1�����t�[���G�ϊ��p�̊֐�
���̃v���O�������g���Ƃ���5�̔z���p�ӂ��Ă���
xr[nx]    : 1�����t�[���G�ϊ��p�̎����f�[�^
xi[nx]    : 1�����t�[���G�ϊ��p�̋����f�[�^
si[nx/2]  : �t�[���G�ϊ��Ŏg���T�C���f�[�^
co[nx/2]  : �t�[���G�ϊ��Ŏg���R�T�C���f�[�^
brv[nx]   : �t�[���G�ϊ��Ŏg���o�^�t���C���Z�f�[�^

nx        : 1�����t�[���G�ϊ��̃f�[�^���i2�̃x�L��j

�֐��F
void FFTInit(int nx, float *si, float *co, unsigned short *brv)
  �t�[���G�ϊ��ɕK�v�ȃT�C���ƃR�T�C���ƃo�^�t���C���Z�f�[�^���쐬����֐�
  FFT�֐����g���O��1�x�g�p���ĕK�v�ȃf�[�^���쐬����

void FFT(int ir, int nx, float *xr, float *xi, float *si, float *co, unsigned short *brv)
  1�����t�[���G�ϊ������s����֐�

void fft2d(int ir, float *fr, float *fi, int nx, int ny)
  2�����t�[���G�ϊ������s����֐�
*/

#include  <stdlib.h>
#include  <math.h>
#define   PI  3.14159265358979

void bitrev(int nx, float *xr, float *xi, unsigned short *brv)
// �o�^�t���C���Z�̓���ւ�
// int   nx;   �f�[�^��
// float *xr;  �����̃f�[�^  xr[nx]
// float *xi;  �����̃f�[�^  xi[nx]
// unsigned short *brv;  �����p�̃f�[�^  brv[nx]
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
// 1�����t�[���G�ϊ�
// int    ir;  ���ϊ�(1)�Ƌt�ϊ�(-1)
// int    nx;  1����FFT�̃f�[�^��
// float *xr;  �����̃f�[�^          xr[nx]
// float *xi;  �����̃f�[�^          xi[nx]
// float *si;  FFT�p�̃T�C���f�[�^   si[nx/2]
// float *co;  FFT�p�̃R�T�C���f�[�^ co[nx/2]
// unsigned short *brv;  FFT�p�̓���ւ��f�[�^ brv[nx]
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
// �����f�[�^�쐬�p�֐�
// int      nx;  �f�[�^��
// unsigned nn;  �����O�̃f�[�^�ԍ�
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
// FFT�p�̃f�[�^�쐬�p�̊֐�
// int    nx;  FFT�̃f�[�^��
// float *si;  �T�C���f�[�^�p�z��   si[nx/2]
// float *co;  �R�T�C���f�[�^�p�z�� co[nx/2]
// unsigned short *brv;  �����f�[�^�p�z�� brv[nx]
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
// 2�����t�[���G�ϊ�
// int    ir;   ���ϊ�(1)�Ƌt�ϊ�(-1)
// float *fr;   2����FFT�̎����̃f�[�^ fr[nx*ny]
// float *fi;   2����FFT�̋����̃f�[�^ fi[nx*ny]
// int    nx;   x�����̃f�[�^��
// int    ny;   y�����̃f�[�^��
{
	int   i, j;
	float *gr, *gi, *si, *co;
	unsigned short *br;

	// x�����̃t�[���G�ϊ�
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

	// y�����̃t�[���G�ϊ�
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
