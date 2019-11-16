/* P2-30qsort.c */

#include  <stdio.h>
#include  <stdlib.h>

int compare(const void *, const void *);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// �N�C�b�N�\�[�g
	qsort(a, n, sizeof(int), compare);
	// ���ʂ̕\��
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

int compare(const void *a, const void *b)
// �N�C�b�N�\�[�g�̔�r�֐�
// const void *a; �\�[�g����1�߂̗v�f
// const void *b; �\�[�g����2�߂̗v�f
// ��r�̂Ƃ��ɁA�L���X�g���Z�q���g���Ĕ�r����v�f�̕ϐ��̌^�ɍ��킹��
{
	if(*(int *)a < *(int *)b)      return -1;
	else if(*(int *)a > *(int *)b) return 1;
	else  return 0;
}
