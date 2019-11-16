/*

x軸の正の方向に検出器を置いている

floatにしてある
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

コンパイルオプションとして
-std=c++11
これを必ずつける（これないとエラーがでる）



*/

#include <stdio.h>
#include <unistd.h>
#include <iostream>
using namespace std;
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
#include "../include/util.h"
#include "../include/util.cpp"
#define DEBUG


typedef struct {
	int detector_size_w;
	int detector_size_h;
	int pinhole_img_w;
	int pinhole_img_h;
	int pinhole_count;
	float img_pixel_size;
	float pinhole_img_pixel_size;
	float rotation_radius;
	float distance_collimator_to_detector;
	float collimator_h;
	float collimator_w;
	float collimator_w_y_axis;
	float collimator_w_z_axis;
	float d_width;
	float d_height;
	float time;
	float photon_num;
	float aperture_degree;
	float aperture_degree_xy;
	float aperture_degree_zx;
} Condition;

void create_fov(std::vector<int> &fov, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num);
void create_pinhole_img(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, int option);
void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, std::vector<int> &fov, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx);
void create_fov(std::vector<int> &fov, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num);
void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num);
void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num);
Eigen::Vector3f rotate_and_move_intersection_edge(Eigen::Vector3f intersection_edge, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx);
void plot_corner(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h, int option);
void connecting_corners_with_lines(std::vector<int> &img, Condition cond, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h);
void fill_rectangle(std::vector<int> &img, Condition cond, int pinhole_num, int width, int height);

// 確認し終わったら消す
void plot_center(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h, int option);


int main()
{
	Condition cond;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.pinhole_img_w = 1024 * 2;
	cond.pinhole_img_h = 512 * 2;
	cond.rotation_radius = 25;
	cond.distance_collimator_to_detector = 7.6;
	cond.collimator_h = 1.;
	// cond.collimator_w = 0.5;
	cond.collimator_w_y_axis = 0.5;
	cond.collimator_w_z_axis = 0.5;
	cond.img_pixel_size = 0.16;
	cond.pinhole_img_pixel_size = 0.04 / 2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.photon_num = 10000000;
	cond.aperture_degree = 24.;
	cond.aperture_degree_xy = 24.;
	cond.aperture_degree_zx = 24.;
	
	std::vector<float> pinhole_theta_xy{ -15., 15., 0., -15., 15. };
	cond.pinhole_count = pinhole_theta_xy.size();
	std::vector<float> pinhole_theta_zx{ 	-9.5f, -9.5f, 0., 9.5f, 9.5f };
	std::vector<float> pinhole_center
	{
	    cond.rotation_radius, - 10.,  3.8,
	    cond.rotation_radius,   10.,  3.8,
	    cond.rotation_radius,    0.,  0.,
	    cond.rotation_radius, - 10., - 3.8,
	    cond.rotation_radius,   10., - 3.8
	};


	// std::vector<float> pinhole_theta_xy{ -24., 24., 0., -24., 24. };
	// cond.pinhole_count = pinhole_theta_xy.size();
	// std::vector<float> pinhole_theta_zx{ 	-9.5f, -9.5f, 0., 9.5f, 9.5f };
	// std::vector<float> pinhole_center
	// {
	//     cond.rotation_radius, - 11.,  4.,
	//     cond.rotation_radius,   11.,  4.,
	//     cond.rotation_radius,    0.,  0.,
	//     cond.rotation_radius, - 11., - 4.,
	//     cond.rotation_radius,   11., - 4.
	// };

	// std::vector<float> pinhole_theta_xy{ 0.};
	// cond.pinhole_count = pinhole_theta_xy.size();
	// std::vector<float> pinhole_theta_zx{ 0.};
	// std::vector<float> pinhole_center
	// {
	//     cond.rotation_radius,    0.,  0.,
	// };

	// std::vector<float> pinhole_theta_xy{ -24.};
	// cond.pinhole_count = pinhole_theta_xy.size();
	// std::vector<float> pinhole_theta_zx{ -9.5};
	// std::vector<float> pinhole_center
	// {
	//     cond.rotation_radius, -3.,  0.,
	// };

	
	/*---- 時間計測開始&条件表示 start ----*/

	/*---- 時間計測開始&条件表示  end -----*/





	// 変数の定義
	std::vector<int> pinhole_img_layer1(cond.pinhole_img_w * cond.pinhole_img_h, -1);
	std::vector<int> pinhole_img_layer3(pinhole_img_layer1.size(), -1);
	std::vector<int> fov(cond.detector_size_w * cond.detector_size_h, -1);
	set_pinhole_img(pinhole_img_layer1, pinhole_img_layer3, fov, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx);

	/*----- 時間計測処理 start-----*/

	/*----- 時間計測処理 end-------*/
}

void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, std::vector<int> &fov, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx)
{
	for(int pinhole_num = 0; pinhole_num < pinhole_theta_xy.size(); pinhole_num++)
	{
		Eigen::Vector3f center;
		center << pinhole_center[coord_num * pinhole_num + X], pinhole_center[coord_num * pinhole_num + Y], pinhole_center[coord_num * pinhole_num + Z];
		float theta_xy = pinhole_theta_xy[pinhole_num];
		float theta_zx = pinhole_theta_zx[pinhole_num];

		create_pinhole_img_layer1(pinhole_img_layer1, cond, center, theta_xy, theta_zx, pinhole_num);
		create_pinhole_img_layer3(pinhole_img_layer3, cond, center, theta_xy, theta_zx, pinhole_num);
		create_fov(fov, cond, center, theta_xy, theta_zx, pinhole_num);
	}

	std::ostringstream ostr;
	{
		ostr << "./result/layer1_" << cond.pinhole_count << "pinhole_square_" << (int)(cond.collimator_w_y_axis * 10) << "-" << (int)(cond.collimator_w_z_axis * 10) << "mm-_int_" << cond.pinhole_img_w << "-" << cond.pinhole_img_h << ".raw";
		std::string write_name = ostr.str(); ostr.str("");
		writeRawFile(write_name.c_str(), pinhole_img_layer1);
	}
	{
		ostr << "./result/layer3_" << cond.pinhole_count << "pinhole_square_" << (int)(cond.collimator_w_y_axis * 10) << "-" << (int)(cond.collimator_w_z_axis * 10) << "mm_int_" << cond.pinhole_img_w << "-" << cond.pinhole_img_h << ".raw";
		std::string write_name = ostr.str(); ostr.str("");
		writeRawFile(write_name.c_str(), pinhole_img_layer3);
	}
	{
		ostr << "./result/fov_" << cond.pinhole_count << "pinhole_square_" << (int)(cond.collimator_w_y_axis * 10) << "-" << (int)(cond.collimator_w_z_axis * 10) << "mm_int_" << cond.detector_size_w << "-" << cond.detector_size_h << ".raw";
		std::string write_name = ostr.str(); ostr.str("");
		writeRawFile(write_name.c_str(), fov);
	}
}

void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 0); }

void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 1); }

void create_fov(std::vector<int> &fov, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num)
{ create_pinhole_img( fov, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 2); }


void create_pinhole_img(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, int option)
{
	int corner_count = 4;
	std::vector<int> corner_i(corner_count * cond.pinhole_count);
	std::vector<int> corner_j(corner_i.size());
	std::vector<int> width{ cond.pinhole_img_w, cond.pinhole_img_w, cond.detector_size_w };
	std::vector<int> height{ cond.pinhole_img_h, cond.pinhole_img_h, cond.detector_size_h };
	std::vector<float> pixel_size_w{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.d_width };
	std::vector<float> pixel_size_h{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.d_height };

	plot_corner(img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, corner_i, corner_j, width[option], height[option], pixel_size_w[option], pixel_size_h[option], option);
	// plot_center(img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, corner_i, corner_j, width[option], height[option], pixel_size_w[option], pixel_size_h[option], option);
	connecting_corners_with_lines(img, cond, pinhole_num, corner_i, corner_j, width[option], height[option], pixel_size_w[option], pixel_size_h[option]);
	fill_rectangle(img, cond, pinhole_num, width[option], height[option]);
	
}


void plot_center(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h, int option)
{
	bool test_bool = pinhole_num == 4;
	int test_plot_num = 2;
	// option { 0: layer1, 1: layer3, 2: fov }

	float collimator_radius_y = cond.collimator_w_y_axis / 2.;
	float collimator_radius_z = cond.collimator_w_z_axis / 2.;
	float aperture_xy = degree_to_rad(cond.aperture_degree_xy);
	float aperture_zx = degree_to_rad(cond.aperture_degree_zx);


	

	std::vector<float> end_x{ static_cast<float>(cond.rotation_radius - cond.collimator_h / 2.),
		 												static_cast<float>(cond.rotation_radius + cond.collimator_h / 2.),
														static_cast<float>(cond.rotation_radius + cond.collimator_h / 2.) };

	Eigen::Vector3f surface_center; surface_center << cond.collimator_h / 2., 0., 0.;
	Eigen::Vector3f surface_center_rot = rotate_and_move_intersection_edge(surface_center, pinhole_center, pinhole_theta_xy, pinhole_theta_zx);

	// ４点の交点のベクトル求める
	Eigen::Vector3f corner_vec;
	corner_vec = surface_center_rot - pinhole_center;

	Eigen::Vector3f on_img;
	std::vector<float> x{ static_cast<float>(cond.rotation_radius - cond.collimator_h / 2.),
		 										static_cast<float>(cond.rotation_radius + cond.collimator_h / 2.),
												static_cast<float>((cond.rotation_radius + cond.distance_collimator_to_detector)) };

	// -------
	if(test_bool)
	{
		cout << "corner_vec << " << corner_vec(X) << ", " << corner_vec(Y) << ", " << corner_vec(Z) << endl;
		float test = tan(corner_vec(Y) / corner_vec(X));
		float tst2 = test * 180 / M_PI;
		cout << "theta: " << tst2 << endl;
	}
	// -------

	
	// ベクトルで画像上の座標を求める
	on_img(X) = x[option];
	float t = (on_img(X) - pinhole_center(X)) / corner_vec(X);
	on_img(Y) = pinhole_center(Y) + t * corner_vec(Y);
	on_img(Z) = pinhole_center(Z) + t * corner_vec(Z);
	int j = transform_img_coordinate_same_axis(on_img(Y), width, pixel_size_w);
	int i = transform_img_coordinate_opposite_axis(on_img(Z), height, pixel_size_h);
	
	// if())
	{
		cout << "on_img << " << on_img(X) << ", " << on_img(Y) << ", " << on_img(Z) << endl;
		cout << "i: " << i << ", j: " << j << endl;
	}
	if(is_out_of_image(i, height, j, width))
	{
		// vector<string> corner_str { "top_left", "top_right", "buttom_right", "buttom_left" };
		// cout << "ERROR : failed to plot corner (pinhole_num: " << pinhole_num << ", corner: " << corner_str[num] << ")" << endl;
		cout << "ERROR" << endl;
		exit(1);
	}
	
	img[width * i + j] = pinhole_num;
}


void plot_corner(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h, int option)
{
	// bool test_bool = pinhole_num == 4;
	bool test_bool = true;
	int test_plot_num = 2;
	// option { 0: layer1, 1: layer3, 2: fov }

	float collimator_radius_y = cond.collimator_w_y_axis / 2.;
	float collimator_radius_z = cond.collimator_w_z_axis / 2.;
	float aperture_xy = degree_to_rad(cond.aperture_degree_xy);
	float aperture_zx = degree_to_rad(cond.aperture_degree_zx);

	int corner_count = 4;
	// 順番：左上，右上，右下，左下
	std::vector<Eigen::Vector3f> pinhole_corner_layer2(corner_count);
	rep(i, corner_count) pinhole_corner_layer2[i](X) = 0.0f;

	pinhole_corner_layer2[0](Y) = - cond.collimator_w_y_axis / 2.;
	pinhole_corner_layer2[1](Y) =   cond.collimator_w_y_axis / 2.;
	pinhole_corner_layer2[2](Y) =   cond.collimator_w_y_axis / 2.;
	pinhole_corner_layer2[3](Y) = - cond.collimator_w_y_axis / 2.;

	pinhole_corner_layer2[0](Z) =  cond.collimator_w_z_axis / 2.;
	pinhole_corner_layer2[1](Z) =  cond.collimator_w_z_axis / 2.;
	pinhole_corner_layer2[2](Z) =  - cond.collimator_w_z_axis / 2.;
	pinhole_corner_layer2[3](Z) =  - cond.collimator_w_z_axis / 2.;

	std::vector<Eigen::Vector3f> pinhole_corner_layer2_rot(pinhole_corner_layer2.size());
	rep(i, corner_count) pinhole_corner_layer2_rot[i] = rotate_and_move_intersection_edge(pinhole_corner_layer2[i], pinhole_center, pinhole_theta_xy, pinhole_theta_zx);

	std::vector<Eigen::Vector3f> pinhole_corner_surface(corner_count);
	std::vector<float> end_x{ static_cast<float>(- cond.collimator_h / 2.),
		 												static_cast<float>(cond.collimator_h / 2.),
														static_cast<float>(cond.collimator_h / 2.) };
	rep(i, corner_count) pinhole_corner_surface[i](X) = end_x[option];
	float edge_y = (cond.collimator_h / 2.) * tan(aperture_xy);
	float edge_z = (cond.collimator_h / 2.) * tan(aperture_zx);

	pinhole_corner_surface[0](Y) = - (cond.collimator_w_y_axis / 2. + edge_y);
	pinhole_corner_surface[1](Y) = cond.collimator_w_y_axis / 2. + edge_y;
	pinhole_corner_surface[2](Y) = cond.collimator_w_y_axis / 2. + edge_y;
	pinhole_corner_surface[3](Y) = - (cond.collimator_w_y_axis / 2. + edge_y);

	pinhole_corner_surface[0](Z) =  cond.collimator_w_z_axis / 2. + edge_z;
	pinhole_corner_surface[1](Z) =  cond.collimator_w_z_axis / 2. + edge_z;
	pinhole_corner_surface[2](Z) =  - (cond.collimator_w_z_axis / 2. + edge_z);
	pinhole_corner_surface[3](Z) =  - (cond.collimator_w_z_axis / 2. + edge_z);

	std::vector<Eigen::Vector3f> pinhole_corner_surface_rot(pinhole_corner_surface.size());
	rep(i, corner_count) pinhole_corner_surface_rot[i] = rotate_and_move_intersection_edge(pinhole_corner_surface[i], pinhole_center, pinhole_theta_xy, pinhole_theta_zx);

	
	// ４点の交点のベクトル求める
	std::vector<Eigen::Vector3f> corner_vec(pinhole_corner_layer2.size());
	rep(i, corner_count) corner_vec[i] = pinhole_corner_surface_rot[i] - pinhole_corner_layer2_rot[i];

	std::vector<Eigen::Vector3f> on_img(pinhole_corner_layer2.size());
	
	std::vector<float> x{ static_cast<float>(cond.rotation_radius - cond.collimator_h / 2.),
		 										static_cast<float>(cond.rotation_radius + cond.collimator_h / 2.),
												static_cast<float>((cond.rotation_radius + cond.distance_collimator_to_detector)) };

	// -------
	if(test_bool)
	{
		cout << "corner_vec << " << corner_vec[test_plot_num](X) << ", " << corner_vec[test_plot_num](Y) << ", " << corner_vec[test_plot_num](Z) << endl;
		float test = tan(corner_vec[test_plot_num](Y) / corner_vec[test_plot_num](X));
		float tst2 = test * 180 / M_PI;
		cout << "theta: " << tst2 << endl;
	}
	// -------

	
	// ベクトルで画像上の座標を求める
	rep(num, corner_count)
	{
		on_img[num](X) = x[option];
		float t = (on_img[num](X) - pinhole_corner_layer2_rot[num](X)) / corner_vec[num](X);
		on_img[num](Y) = pinhole_corner_layer2_rot[num](Y) + t * corner_vec[num](Y);
		on_img[num](Z) = pinhole_corner_layer2_rot[num](Z) + t * corner_vec[num](Z);
		int j = corner_j[pinhole_num * corner_count + num] = transform_img_coordinate_same_axis(on_img[num](Y), width, pixel_size_w);
		int i = corner_i[pinhole_num * corner_count + num] = transform_img_coordinate_opposite_axis(on_img[num](Z), height, pixel_size_h);
		
		if(test_bool && num == test_plot_num)
		{
			cout << "on_img << " << on_img[test_plot_num](X) << ", " << on_img[test_plot_num](Y) << ", " << on_img[test_plot_num](Z) << endl;
			cout << "i: " << i << ", j: " << j << endl;
		}
		if(is_out_of_image(i, height, j, width))
		{
			vector<string> corner_str { "top_left", "top_right", "buttom_right", "buttom_left" };
			cout << "ERROR : failed to plot corner (pinhole_num: " << pinhole_num << ", corner: " << corner_str[num] << ")" << endl;
			exit(1);
		}
		
		img[width * i + j] = pinhole_num;
	}
}

	
void connecting_corners_with_lines(std::vector<int> &img, Condition cond, int pinhole_num, std::vector<int> &corner_i, std::vector<int> &corner_j, int width, int height, float pixel_size_w, float pixel_size_h)
{
	
	// 4かいループ
	// 左上->右上、右上->右下、右下->左下、左下->左上
	int num_of_sides = 4;
	int corner_count = 4;
	rep(num, num_of_sides)
	{
		Eigen::Vector2f start, end, corner_vec, sp;
		start(X) = transform_origin_coordinate_same_axis(corner_j[pinhole_num * corner_count + num], width, pixel_size_w);
		start(Y) = transform_origin_coordinate_opposite_axis(corner_i[pinhole_num * corner_count + num], height, pixel_size_h);

		int index = (num + 1 == corner_count) ? pinhole_num * corner_count : pinhole_num * corner_count + num + 1;
		// int index = num + 1;
		end(X) = transform_origin_coordinate_same_axis(corner_j[index], width, pixel_size_w);
		end(Y) = transform_origin_coordinate_opposite_axis(corner_i[index], height, pixel_size_h);

		// ベクトル作る
		corner_vec = calculate_unit_vector_2f(start, end);

		sp = start;
		sp(X) /= pixel_size_w;
		sp(Y) /= pixel_size_h;

		// その始点から終点まで点を売っていく（間隔：ピクセルサイズ）
		int counter = 0;
		int before_i, before_j;
		for(int sp_count = 0; counter < 1; sp_count++)
		{

			int j = transform_img_coordinate_same_axis(sp(X), width);
			int i = transform_img_coordinate_opposite_axis(sp(Y), height);
			if(i == before_i && j == before_j) { sp += corner_vec; continue; }
			
			// 終点まで達しば場合は終了（多少のズレがあるため、大まかは判断ししている）
			bool recognize_value;
			if(num == 0 || num == 2) recognize_value = img[(i - 1) * width + j] == pinhole_num || img[(i + 1) * width + j] == pinhole_num || img[i  * width + j] == pinhole_num;
			else recognize_value = img[i * width + j - 1] == pinhole_num || img[i  * width + j + 1] == pinhole_num || img[i  * width + j] == pinhole_num;

			bool around_end_point = ((i - 1 == corner_i[index] && j - 1 == corner_j[index]) ||
															 (i - 1 == corner_i[index] && j == corner_j[index]) ||
															 (i - 1 == corner_i[index] && j + 1 == corner_j[index]) ||
															 (i == corner_i[index] && j - 1 == corner_j[index]) ||
															 (i == corner_i[index] && j == corner_j[index]) ||
															 (i == corner_i[index] && j + 1 == corner_j[index]) ||
															 (i + 1 == corner_i[index] && j - 1 == corner_j[index]) ||
															 (i + 1 == corner_i[index] && j == corner_j[index]) ||
															 (i + 1 == corner_i[index] && j + 1 == corner_j[index]));
			
			bool reach_end_point = sp_count >= 5 && recognize_value && around_end_point;
			if(reach_end_point) counter++;
			
			img[i * width + j] = pinhole_num;
			sp += corner_vec;
			before_i = i;
			before_j = j;
			sp_count++;
		}
	}
}


void fill_rectangle(std::vector<int> &img, Condition cond, int pinhole_num, int width, int height)
{
	std::vector<int> pinhole_exist_i(height, 0);
	rep(i, height)
	{
		// 穴を埋める部分があるかどうかを確認するために、-1とpinhole_numが交互にあるかを確認
		rep(j, width)
		{
			if(pinhole_exist_i[i] == 1 || pinhole_exist_i[i] == 3) { if(img[i * width + j] == pinhole_num) pinhole_exist_i[i]++; }
			else if(pinhole_exist_i[i] == 0 || pinhole_exist_i[i] == 2 || pinhole_exist_i[i] == 4) { if(img[i * width + j] == -1) pinhole_exist_i[i]++; }
		}
		
	}

	rep(i, height)
	{
		if(pinhole_exist_i[i] <= 4) continue;

		int corner_count = 0;
		rep(j, width)
		{
			if(corner_count == 1 || corner_count == 3) { if(img[i * width + j] == pinhole_num) corner_count++; }
			else if(corner_count == 0 || corner_count == 2 || corner_count == 4) { if(img[i * width + j] == -1) corner_count++; }

			if(corner_count == 3) img[i * width + j] = pinhole_num;
		}
	}
}

Eigen::Vector3f rotate_and_move_intersection_edge(Eigen::Vector3f intersection_edge, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx)
{
	Eigen::Matrix3f rot_xy, rot_zx;
	Eigen::Vector3f axis; axis << 0, 0, 1;
	rot_xy = Eigen::AngleAxisf( degree_to_rad(pinhole_theta_xy), axis );

	axis << 0, 1, 0;
	rot_zx = Eigen::AngleAxisf( degree_to_rad(pinhole_theta_zx), axis );

	return rot_zx * rot_xy * intersection_edge + pinhole_center;
}
