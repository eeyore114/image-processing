/* P2-30qsort.c */

#include  <stdio.h>
#include  <stdlib.h>

int compare(const void *, const void *);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// クイックソート
	qsort(a, n, sizeof(int), compare);
	// 結果の表示
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

int compare(const void *a, const void *b)
// クイックソートの比較関数
// const void *a; ソートする1つめの要素
// const void *b; ソートする2つめの要素
// 比較のときに、キャスト演算子を使って比較する要素の変数の型に合わせる
{
	if(*(int *)a < *(int *)b)      return -1;
	else if(*(int *)a > *(int *)b) return 1;
	else  return 0;
}
