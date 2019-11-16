/* P2-31rand.c */

#include  <stdio.h>
#include  <stdlib.h>

main( )
{
	int    i, n = 5;
	double max = RAND_MAX;

	// 疑似乱数ジェネレータの初期化
	srand(12345);

	// 疑似乱数の表示（整数値）
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %d \n", i+1, rand() );

	// 0から99の間の乱数（整数値）
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %d \n", i+1, (int)(100.*rand()/max) );

	// 0以上1未満の間の乱数（実数値）
	for(i = 0 ; i < n ; i++)
		printf("[%2d] %f \n", i+1, (double)rand()/max );
}
