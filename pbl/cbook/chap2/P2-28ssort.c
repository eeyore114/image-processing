/* P2-28ssort.c */

#include  <stdio.h>

void ssort(int *, int);

main( )
{
	int i, n = 5;
	int a[5] = {3, 5, 1, 4, 2};
	// 選択ソート
	ssort(a, n);
	// 結果の表示
	for(i = 0 ; i < n ; i++)
		printf("a[%d]=%d \n", i, a[i]);
}

void  ssort(int *a, int n)
// 選択ソートの関数
// int  *a; ソートする配列
// int  n;  配列の要素数
{
	int  i, j, buff;
	for(i = 0 ; i < n-1 ; i++) {
		for(j = i+1 ; j < n ; j++) {
			if(a[i] > a[j]) { // 交換（スワップ）
				buff = a[i];
				a[i] = a[j];
				a[j] = buff;
			}
		}
	}
}
