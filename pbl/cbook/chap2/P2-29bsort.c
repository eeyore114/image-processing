/* P2-29bsort.c */

#include  <stdio.h>

void bsort(int *, int);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// �o�u���\�[�g
	bsort(a, n);
	// ���ʂ̕\��
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

void  bsort(int *a, int n)
// �o�u���\�[�g�̊֐�
// int  *a; �\�[�g����z��
// int  n;  �z��̗v�f��
{
	int  i, j, buff;
	for(i = 0 ; i < n-1 ; i++) {
		for(j = n-2 ; j >= i ; j--) {
			if(a[j] > a[j+1]) { // �����i�X���b�v�j
				buff = a[j];
				a[j] = a[j+1];
				a[j+1] = buff;
			}
		}
	}
}
