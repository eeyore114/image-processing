/*  P5-03backproj.c  */

#include  <math.h>
#define   PI  3.14159265358979

void BackProjection(int pp, float *img, int nx, int ny, double plx, double ply, float *prj, int px, int pa, double pl)
// �t���e���s���֐�
// int    pp;   �t���e��PI�ōs����2*PI�ōs�����i1 or 2�j
// float  *img; �č\�������摜�f�[�^
// int    nx;   �摜�̃}�g���N�X�T�C�Y�ix�����j
// int    ny;   �摜�̃}�g���N�X�T�C�Y�iy�����j
// double plx;  �摜�̃s�N�Z�������ix�����Fcm�j
// double ply;  �摜�̃s�N�Z�������iy�����Fcm�j
// float  *prj; ���e�f�[�^
// int    px;   ���e�f�[�^�̓��a�����̃f�[�^��
// int    pa;   ���e�f�[�^�̊p�x�����̃f�[�^��
// double pl    ���e�f�[�^�̓��a�����̃s�N�Z�������icm�j
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
