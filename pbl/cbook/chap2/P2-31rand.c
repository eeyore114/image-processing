/* P2-31rand.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	int    i, n = 5;
	double max = RAND_MAX;

	// �^�������W�F�l���[�^�̏�����
	srand(12345);

	// �^�������̕\���i�����l�j
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %d \n", i+1, rand() );

	// 0����99�̊Ԃ̗����i�����l�j
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %d \n", i+1, (int)(100.*rand()/max) );

	// 0�ȏ�1�����̊Ԃ̗����i�����l�j
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %f \n", i+1, (double)rand()/max );
}
