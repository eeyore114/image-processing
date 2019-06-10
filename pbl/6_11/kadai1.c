#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//課題1:x-y座標平面上の点を表す構造体
 typedef struct
 {
	double x;
	double y;
	
} coordinate_xy;

//課題2
void show_xy(coordinate_xy coordinate);
//課題3 : 点P,Qを受け取り、x座標、y座標をそれぞれ加算して値を返す。
coordinate_xy add_xy(coordinate_xy P,coordinate_xy Q);
//課題4
int Inner_Product(coordinate_xy P, coordinate_xy Q);
//課題5 : 二直線のなす角を度数法で返す
int calculate_degree(coordinate_xy P,coordinate_xy Q,int inner_product);
//課題6 : 座標平面上の点を原点を基準にn倍拡大
void multi_n_times(coordinate_xy *P, double n);

int main()
{
	coordinate_xy P;
	P.x = 4.0f;
	P.y = 1.0f;
	coordinate_xy Q;
	Q.x = 0.0f;
	Q.y = 5.0f;
	
	//課題2
	printf("task2\n");
	show_xy(P);
	printf("-----------------\n");

	//課題3
	printf("task3\n");
	coordinate_xy add = add_xy(P, Q);
	show_xy(add);
	printf("-----------------\n");

	//課題4
	printf("task4\n");
	int inner_product = Inner_Product(P, Q);
	printf("inner_product = %d\n", inner_product);
	printf("-----------------\n");	
	
	//課題5
	printf("task5\n");
	int theta = calculate_degree(P, Q, inner_product);
	printf("theta = %d\n", theta);
	printf("-----------------\n");	
	
	//課題6
	printf("task6\n");
	double n = 0.5f;
	multi_n_times(&P, n);
	show_xy(P);
	printf("-----------------\n");

	//課題7
	printf("task7\n");
	
	printf("-----------------\n");	


}

//課題2
void show_xy(coordinate_xy coordinate)
{
	printf("x = %.1f\n", coordinate.x);
	printf("y = %.1f\n", coordinate.y);
}

//課題3:点P,Qを受け取り、x座標、y座標をそれぞれ加算して値を返す。
coordinate_xy add_xy(coordinate_xy P, coordinate_xy Q)
{	
	coordinate_xy add;
	add.x = P.x + Q.x;
	add.y = P.y + Q.y;
	return add;
}


//課題4
int Inner_Product(coordinate_xy P, coordinate_xy Q)
{
	int inner_product = P.x * Q.x + P.y * Q.y;
	return inner_product;
}


//課題5 : 二直線のなす角を度数法で返す
int calculate_degree(coordinate_xy P,coordinate_xy Q, int inner_product)
{
	double P_distance 	  = sqrt(pow(P.x, 2.0) + pow(P.y, 2.0));
	double Q_distance     = sqrt(pow(Q.x, 2.0) + pow(Q.y, 2.0));
	double cos_theta_rad  = inner_product / (P_distance * Q_distance);
	double theta_rad 	  = acos(cos_theta_rad);
	//弧度法を度数法に変換
	int    theta_degree   = (int)(theta_rad * 180.0f / M_PI);
	return theta_degree;
}


//課題6 : 座標平面上の点を原点を基準にn倍拡大
void multi_n_times(coordinate_xy *P, double n)
{
	P->x *= n;
	P->y *= n;
}


