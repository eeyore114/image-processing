/* P2-29bsort.c */

#include  <stdio.h>

void bsort(int *, int);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// バブルソート
	bsort(a, n);
	// 結果の表示
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

void  bsort(int *a, int n)
// バブルソートの関数
// int  *a; ソートする配列
// int  n;  配列の要素数
{
	int  i, j, buff;
	for(i = 0 ; i < n-1 ; i++) {
		for(j = n-2 ; j >= i ; j--) {
			if(a[j] > a[j+1]) { // 交換（スワップ）
				buff = a[j];
				a[j] = a[j+1];
				a[j+1] = buff;
			}
		}
	}
}
