/*

2019/7/28　作成
2直線の交点を求める関数

ナイフエッジの交点を求めるために作成したため，3次元のEigenで定義しているが，計算はxy平面のみで実装している．

*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
// #define DEBUG

typedef struct {
	// y = (grad)x + intercept とする
	float grad;
	float intercept;
} Equation1d;

Equation1d calc_liner_equation(Eigen::Vector3f start, Eigen::Vector3f end)
{
	Equation1d equation;
	equation.grad = (end(1) - start(1)) / (end(0) - start(0));
	equation.intercept = start(1) - equation.grad * start(0);
	return equation;
}

Eigen::Vector3f calc_cross_point(Eigen::Vector3f left_start, Eigen::Vector3f left_end, Eigen::Vector3f right_start, Eigen::Vector3f right_end)
{
	Equation1d left, right;
	left = calc_liner_equation(left_start, left_end);
	right = calc_liner_equation(right_start, right_end);
	Eigen::Vector3f cross_point;
	// zはそのまま，xy平面で2直線の交点を求める
	cross_point(2) = left_start(2);
	cross_point(0) = (right.intercept - left.intercept) / (right.grad - left.grad);
	cross_point(1) = right.grad * cross_point(0) + right.intercept;
	return cross_point;
}

void test()
{

	Eigen::Vector3f left_start, left_end, right_start, right_end;
	float z = 5.;
	left_start << -2., -1, z;
	left_end << -1, 0, z;
	right_start << 2, -3, z;
	right_end << 1, 0, z;

	Eigen::Vector3f cross_point = calc_cross_point(left_start, left_end, right_start, right_end);
	std::cout << "cross_point" << std::endl;
	std::cout << cross_point << std::endl;
}

int main()
{
	test();
}
