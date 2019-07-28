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


typedef struct {
	// y = (grad)x + intercept とする
	float grad;
	float intercept;
} Equation1d;


void create_collimator_img(std::vector<float> &detector, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz);
void set_cross_point(std::vector<Eigen::Vector3f> &cross_point, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz);
Equation1d calc_liner_equation(Eigen::Vector3f start, Eigen::Vector3f end);
Eigen::Vector3f calc_cross_point(Eigen::Vector3f left_start, Eigen::Vector3f left_end, Eigen::Vector3f right_start, Eigen::Vector3f right_end);




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
	create_collimator_img(detector, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_yz);

	// 関数呼び出し
	// launch_test_gradient_pinhole(detector, pinhole_theta_xy, pinhole_theta_yz, pinhole_center, cond);

	/*----- 時間計測処理 start-----*/

	/*----- 時間計測処理 end-------*/
}

void create_collimator_img(std::vector<float> &detector, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz)
{
	/*
		- 2つのエッジの交点を求める
		- ピンホールの表面上のピクセルをひとつづつ見て，内積を比較
			- 表面の座標を見て，距離が近い方と内積を比較　
	*/

	/* --- Eigen の配列を作れるらしい --- */
	// std::vector<Eigen::Vector3f> v(2);
	// v[0] << 1, 2, 3;
	// v[1] << 4, 5, 6;
	// // for(int i = 0; i < 2; i++) std::cout << v[i] << std::endl;
	//
	// std::cout << v[0](1) << '\n';
	/* --- Eigen の配列を作れるらしい --- */








	std::vector<Eigen::Vector3f> cross_point_edge(pinhole_theta_xy.size());
	set_cross_point(cross_point_edge, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_yz);









	/*---- ナイフエッジの交点を求める -----*/
	// Eigen::Vector3f cross_point_edge;
	// {
	// 	float x1 = cond.collimator_w / 2.;
	// 	float y1 = -cond.collimator_h / 2.;
	// 	float x2 = cond.collimator_w / 2. + (cond.collimator_h / 2.) * tan(M_PI / 6.);
	// 	float y2 = -cond.collimator_h;
	// 	float t = - x1 / (x2 - x1);
	// 	float y = y1 + t * (y2 - y1);
	//
	// 	cross_point_edge << 0, y, 0;
	// }




	/*
	・決まった座標がないから今までみたいに求めることができない．
	・tの求め方がわからない


	*/





}


void set_cross_point(std::vector<Eigen::Vector3f> &cross_point, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz)
{
	/*
	pinholeの数だけループ
	4点を求める
	func : calc_cross_point

	*/
	for(int i = 0; i < pinhole_theta_xy.size(); i++)
	{
		Eigen::Vector3f left_start, left_end, right_start, right_end;
		// 4点を求める

		cross_point[i] = calc_cross_point(left_start, left_end, right_start, right_end);
	}
}


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
