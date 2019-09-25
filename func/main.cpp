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
#include "head.h"
#include "head.cpp"
// #define DEBUG

int main()
{
	std::cout << coord_num << '\n';
	float theta = rad_to_degree(M_PI / 2.0f);
	std::cout << "theta = " << theta << '\n';
}
