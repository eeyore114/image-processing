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
#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/Geometry"
#include "../include/fileio.h"
#include "../include/Mersenne_twister.h"
#include "../include/Mersenne_twister.cpp"
// #define DEBUG

enum Coordinate { X, Y, Z, coord_num };

typedef struct {
	int detector_size_w;
	int detector_size_h;
	float img_pixel_size;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float d_width;
	float d_height;
	float time;
	float photon_num;
	float knife_edge_theta_degree;
} Condition;

int main()
{
	Condition cond;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 13;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 5;
	cond.img_pixel_size = 0.2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.photon_num = 10000000;
	cond.knife_edge_theta_degree = 24;
	std::vector<float> pinhole_theta_xy{ -24., -15., 15., 24., -9., 0., 9., -24., -15., 15., 24. };
	std::vector<float> pinhole_theta_yz{ 	 9.,   9.,  9.,  9.,  0., 0., 0., - 9., - 9., -9., -9. };
	float pc_y = - (cond.rotation_radius + cond.height_collimator / 2.);
	std::vector<float> pinhole_center
	{ -10., pc_y,  4.2f,
		- 3., pc_y,  4.2f,
		  3., pc_y,  4.2f,
		 10., pc_y,  4.2f,
		-6.2f, pc_y,  0.,
		  0., pc_y,  0.,
		 6.2f, pc_y,  0.,
		-10., pc_y, -4.2f,
 		- 3., pc_y, -4.2f,
 		  3., pc_y, -4.2f,
 		 10., pc_y, -4.2f,
										};

	/*---- 時間計測開始&条件表示 start ----*/

	/*---- 時間計測開始&条件表示  end -----*/

	// 変数の定義
	std::vector<float> detector(cond.detector_size_w * cond.detector_size_h, 0.);

	// 関数呼び出し
	// launch_test_gradient_pinhole(detector, pinhole_theta_xy, pinhole_theta_yz, pinhole_center, cond);

	/*----- 時間計測処理 start-----*/

	/*----- 時間計測処理 end-------*/
}

void create_collimator_img()
{
	/*
		- 2つのエッジの交点を求める
		- ピンホールの表面上のピクセルをひとつづつ見て，内積を比較
			- 表面の座標を見て，距離が近い方と内積を比較　
	*/


	
}
