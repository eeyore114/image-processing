/*

floatにしてある
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

コンパイルオプションとして
-std=c++11
これを必ずつける（これないとエラーがでる）
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
#include "fileio.h"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define DEBUG

int main()
{
	int num = 2;
	Eigen::Vector3f v;
	v << 1, 2, 3;
	std::cout << v(num) << '\n';
}
