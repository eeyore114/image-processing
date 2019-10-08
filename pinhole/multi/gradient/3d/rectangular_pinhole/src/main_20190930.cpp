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
bool pass_rectangular_pinhole(Eigen::Vector2f intersection_to_layer, Eigen::Vector2f base, float aperture, Condition cond);
Eigen::Vector3f rotate_and_move_intersection_edge(Eigen::Vector3f intersection_edge, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx);
float degree_to_rad(float theta_degree);
float rad_to_degree(float theta);


bool pass_rectangular_pinhole_y(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond);
bool pass_rectangular_pinhole_z(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond);


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
	// cond.pinhole_img_pixel_size = 0.08 / 4;
	cond.pinhole_img_pixel_size = 0.04 / 2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.photon_num = 10000000;
	cond.aperture_degree_xy = 24.;
	cond.aperture_degree_zx = 24.;
	std::vector<float> pinhole_theta_xy{ -24., 24., 0., -24., 24. };
	// std::vector<float> pinhole_theta_xy(11, 0.);
	std::vector<float> pinhole_theta_zx{ 	-9.5f, -9.5f, 0., 9.5f, 9.5f };
	// std::vector<float> pinhole_theta_zx(11, 0.);
	std::vector<float> pinhole_center
	{
	    cond.rotation_radius, - 11.,  4.,
	    cond.rotation_radius,   11.,  4.,
	    cond.rotation_radius,    0.,  0.,
	    cond.rotation_radius, - 11., - 4.,
	    cond.rotation_radius,   11., - 4.
	};

	// {
	//     cond.rotation_radius, - 12.,  4.,
	//     cond.rotation_radius, -  3, 4.7f,
	//     cond.rotation_radius,    3., 4.6f,
	//     cond.rotation_radius,   13.2f,  4.,
	//     cond.rotation_radius, - 6.5,  0.,
	//     cond.rotation_radius,    0.,  0.,
	//     cond.rotation_radius,   6.7,  0.,
	//     cond.rotation_radius, - 11., - 4.,
	//     cond.rotation_radius, -  3., - 4.5,
	//     cond.rotation_radius,    3., - 4.5,
	//     cond.rotation_radius,   11., - 4.
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
	/*
	実装のアルゴリズム
	https://paper.dropbox.com/doc/--AjQNF888gr4SxY2nmhuvsAERAQ-zfZOaNRJYJ2QN2vVwrZTU
	*/

	// option { 0: layer1, 1: layer3, 2: fov }
	std::vector<int> width{ cond.pinhole_img_w, cond.pinhole_img_w, cond.detector_size_w };
	std::vector<int> height{ cond.pinhole_img_h, cond.pinhole_img_h, cond.detector_size_h };
	std::vector<float> pixel_size_w{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.d_width };
	std::vector<float> pixel_size_h{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.d_height };
	std::vector<float> sign{ 1., -1., - 1. };

	float collimator_radius_y = cond.collimator_w_y_axis / 2.;
	float collimator_radius_z = cond.collimator_w_z_axis / 2.;
	float aperture_xy = degree_to_rad(cond.aperture_degree_xy);
	float aperture_zx = degree_to_rad(cond.aperture_degree_zx);
	Eigen::Vector3f intersection_edge_y;  intersection_edge_y << sign[option] * collimator_radius_y / tan(aperture_xy), 0., 0.;
	Eigen::Vector3f intersection_edge_z;  intersection_edge_z << sign[option] * collimator_radius_z / tan(aperture_zx), 0., 0.;

	Eigen::Vector3f intersection_edge_y_rot,intersection_edge_z_rot ;
	intersection_edge_y_rot = rotate_and_move_intersection_edge(intersection_edge_y, pinhole_center, pinhole_theta_xy, pinhole_theta_zx);
	intersection_edge_z_rot = rotate_and_move_intersection_edge(intersection_edge_z, pinhole_center, pinhole_theta_xy, pinhole_theta_zx);


	std::vector<float> x{ static_cast<float>(cond.rotation_radius - cond.collimator_h / 2.),
		 										static_cast<float>(cond.rotation_radius + cond.collimator_h / 2.),
												static_cast<float>((cond.rotation_radius + cond.distance_collimator_to_detector)) };
	for (int i = 0; i < height[option]; i++)
	{
		for (int j = 0; j < width[option]; j++)
		{
			float y = (j - (width[option] - 1.) / 2.) * pixel_size_w[option];
			float z = ((height[option] - 1.) / 2. - i) * pixel_size_h[option];
			Eigen::Vector3f layer; layer << x[option], y, z;

			Eigen::Vector3f base_y; base_y = pinhole_center - intersection_edge_y_rot;
			Eigen::Vector3f intersection_to_layer_y; intersection_to_layer_y = layer - intersection_edge_y_rot;
			bool is_in_focus_y = pass_rectangular_pinhole_y(intersection_to_layer_y, base_y, aperture_xy, cond);
			if ( !is_in_focus_y ) continue;

			Eigen::Vector3f base_z; base_z = pinhole_center - intersection_edge_z_rot;
			Eigen::Vector3f intersection_to_layer_z; intersection_to_layer_z = layer - intersection_edge_z_rot;
			bool is_in_focus_z = pass_rectangular_pinhole_z(intersection_to_layer_z, base_z, aperture_zx, cond);
			if ( !is_in_focus_z ) continue;

			img[width[option] * i + j] = pinhole_num;
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

bool pass_rectangular_pinhole_y(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond)
{
	Eigen::Vector2f intersection_to_layer_xy; intersection_to_layer_xy << intersection_to_layer(0), intersection_to_layer(1);
	Eigen::Vector2f base_xy; base_xy << base(0), base(1);
	return pass_rectangular_pinhole(intersection_to_layer_xy, base_xy, aperture, cond);
}


bool pass_rectangular_pinhole_z(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond)
{
	Eigen::Vector2f intersection_to_layer_zx; intersection_to_layer_zx << intersection_to_layer(0), intersection_to_layer(2);
	Eigen::Vector2f base_zx; base_zx << base(0), base(2);
	return pass_rectangular_pinhole(intersection_to_layer_zx, base_zx, aperture, cond);
}

bool pass_rectangular_pinhole(Eigen::Vector2f intersection_to_layer, Eigen::Vector2f base, float aperture, Condition cond)
{
	float dot = intersection_to_layer.dot(base);
	float norm1 = intersection_to_layer.norm();
	float norm2 = base.norm();
	if ( acosf( dot / (norm1 * norm2) ) > aperture ) return false;

	return true;
}


float degree_to_rad(float theta_degree) {	return theta_degree * M_PI / 180.; }

float rad_to_degree(float theta) { return theta * 180. / M_PI; }
