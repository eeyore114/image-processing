/* P2-28ssort.c */

#include  <stdio.h>

void ssort(int *, int);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// �I���\�[�g
	ssort(a, n);
	// ���ʂ̕\��
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

void  ssort(int *a, int n)
// �I���\�[�g�̊֐�
// int  *a; �\�[�g����z��
// int  n;  �z��̗v�f��
{
	int  i, j, buff;
	for(i = 0 ; i < n-1 ; i++) {
		for(j = i+1 ; j < n ; j++) {
			if(a[i] > a[j]) { // �����i�X���b�v�j
				buff = a[i];
				a[i] = a[j];
				a[j] = buff;
			}
		}
	}
}
