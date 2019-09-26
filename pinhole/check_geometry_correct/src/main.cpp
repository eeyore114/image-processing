/*

とりあえず矩形でやってみる

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

コンパイルオプションとして
-std=c++11
これを必ずつける（これないとエラーがでる）

通った場所の特定をしようとしたらうまく行かない（origin_fov）

交点逆かも
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
#include "../include/util.h"
#include "../include/util.cpp"
#include "../include/Mersenne_twister.h"
#include "../include/Mersenne_twister.cpp"
// #define DEBUG

enum PinholeShape { CIRCLE, RECTANGLE };

typedef struct {
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int surface_source_w;
	int surface_source_h;
	int img_w;
	int img_h;
	int img_d;
	int update_count;
	int photon_scale;
	int photon_num_efficiency_map;
	int max_scattering_csount;
	int pinhole_count;
	int pinhole_img_w;
	int pinhole_img_h;
	int pinhole_shape;
	int thread_num;
	int ray_num;
	int max_scattering_count;
	bool is_voxel;
	bool has_efficiency_map;
	float cut_off_energy;
	float rotation_radius;
	float distance_collimator_to_detector;
	float collimator_h;
	/*------------------*/
	// 円形ピンホールの場合
	float collimator_w;
	float aperture_degree;
	// 矩形ピンホールの場合
	float collimator_w_y_axis;
	float collimator_w_z_axis;
	float aperture_degree_xy;
	float aperture_degree_zx;
	/*------------------*/
	float img_pixel_size;
	float detector_pixel_size_w;
	float detector_pixel_size_h;
	float pinhole_img_pixel_size;
	float time;
} Condition;

typedef struct {
	std::string img_file;
	std::string result_file;
	std::string log_file;
	std::string date_directory;
	std::string medium;
	std::string init_position;
	std::string efficiency_map_name;
	std::string read_activity_map_name;
	std::string read_absorp_map_name;
	std::string read_surface_source_name;
	std::string read_knife_edge_layer1_name;
	std::string read_knife_edge_layer3_name;
	std::string write_reconstruct_img_name;
	std::string write_detector_name;
	std::string write_primary_detector_name;
	std::string write_efficiency_filter_name;
	std::string write_efficiency_correction_name;
	std::string write_efficiency_map_standardization_name;
	std::string read_fov_name;

	std::string write_origin_fov_name;

	// 感度補正から行う場合
	std::string read_detector_name;
	std::string read_primary_detector_name;
} HostCondition;

std::string makeLogDirectory();
std::string showCurrentTime();
void outputLogInit(Condition cond, HostCondition cond_host);
void outputLogLast(Condition cond, HostCondition cond_host);
void launchCheck(std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float> pinhole_center, Condition cond, HostCondition cond_host);
void launchCheckGeometryCorrect(std::vector<int>& fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, Condition cond, HostCondition cond_host);
void output_test_result(bool is_correct);
bool CheckGometryCorrect(std::vector<float>& img, std::vector<int>& fov, std::vector<int>& origin_fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, Condition cond);
void launchCreateGeometry(std::vector<int>& fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float> pinhole_center, Condition cond, HostCondition cond_host);

// 画像作成
void create_fov(std::vector<int> &fov, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z);
void create_pinhole_img(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, int option, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z);
void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, std::vector<int> &fov, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z);
void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z);
void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z);
bool pass_rectangular_pinhole(Eigen::Vector2f intersection_to_layer, Eigen::Vector2f base, float aperture, Condition cond);
Eigen::Vector3f rotate_and_move_intersection_edge(Eigen::Vector3f intersection_edge, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx);
bool pass_rectangular_pinhole_y(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond);
bool pass_rectangular_pinhole_z(Eigen::Vector3f intersection_to_layer, Eigen::Vector3f base, float aperture, Condition cond);
void showOriginFOV(std::vector<int> &origin_fov, Eigen::Vector3f cross_point, Eigen::Vector3f on_detector, Eigen::Vector3f d, Condition cond, int pinhole_num, int axis);


int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.img_d = 128;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.surface_source_w = 2048;
	cond.surface_source_h = 1024;
	cond.rotation_radius = 25;
	cond.photon_scale = 30000;
	cond.photon_num_efficiency_map = 1E6;
	cond.pinhole_img_w = 2048;
	cond.pinhole_img_h = 1024;
	cond.thread_num = 4;
	cond.ray_num = 7;
	cond.max_scattering_count = 5;
	cond.cut_off_energy = 30.;
	cond.distance_collimator_to_detector = 7.6;
	cond.collimator_h = 1.;
	/*-----------------------------*/
	// 円形ピンホールの場合
	// cond.collimator_w = 0.5;
	// cond.aperture_degree = 24.0f;
	// 矩形ピンホールの場合
	cond.collimator_w_y_axis = 0.2;
	cond.collimator_w_z_axis = 0.5;
	cond.aperture_degree_xy = 20.;
	cond.aperture_degree_zx = 30.;
	/*-----------------------------*/
	cond.pinhole_shape = RECTANGLE;
	cond.update_count = 1;
	cond.img_pixel_size = 0.15;
	cond.detector_pixel_size_w = 0.08;
	cond.detector_pixel_size_h = 0.08;
	cond.pinhole_img_pixel_size = 0.02;
	cond.has_efficiency_map = true;
	// std::vector<float> pinhole_theta_xy{ -24., 24., 0., -24., 24. };
	// std::vector<float> pinhole_theta_zx{ -9.5f, -9.5f, 0., 9.5f, 9.5f };
	// cond.pinhole_count = pinhole_theta_xy.size();
	//
	// std::vector<float> pinhole_center
	// {
	//     cond.rotation_radius, - 11.,  4.,
	//     cond.rotation_radius,   11.,  4.,
	//     cond.rotation_radius,    0.,  0.,
	//     cond.rotation_radius, - 11., - 4.,
	//     cond.rotation_radius,   11., - 4.
	// };


	std::vector<float> pinhole_theta_xy{ -25.6, -14.8, 0., 13., 25.6, -25.6, -14.8, 0., 13., 25.6 };
	std::vector<float> pinhole_theta_zx{ -7.9f, -7.9f, -7.9f,  -7.9f, -7.9f, 7.9f, 7.9f, 7.9f, 7.9f, 7.9f };
	cond.pinhole_count = pinhole_theta_xy.size();
	std::vector<float> pinhole_center
	{
	    cond.rotation_radius, - 12,  3.5,
			cond.rotation_radius,    -5.5,  3.5,
			cond.rotation_radius,    0.,  3.5,
	    cond.rotation_radius,    5.5,  3.5,
			cond.rotation_radius,   12,  3.5,
	    cond.rotation_radius, - 12, - 3.5,
			cond.rotation_radius,    -5.5, - 3.5,
			cond.rotation_radius,     0., - 3.5,
	    cond.rotation_radius,    5.5, - 3.5,
	    cond.rotation_radius,   12, - 3.5
	};

	HostCondition cond_host;
	cond_host.init_position = "voxel"; //"voxel" or "origin"
	cond.is_voxel = cond_host.init_position == "voxel";

	cond_host.medium = "Shepp"; //"Shepp", "Brain", "Sphere" など（大文字始まり）
	cond_host.img_file = "./img/";
	cond_host.result_file = "./result/";
	cond_host.log_file = "./log/";
	cond_host.date_directory = makeLogDirectory();
	cond_host.efficiency_map_name = "efficiency_map_rectangle_float_512-256.raw";
	// cond_host.read_activity_map_name = "brain_float_128-128-128.raw";
	cond_host.read_activity_map_name = "Shepp_float_128-128-128.raw";
	cond_host.read_absorp_map_name = "Shepp_absorp_map_float_64-64-64.raw";
	cond_host.read_surface_source_name = "surfaceSource_int_2048-1024.raw";
	cond_host.write_reconstruct_img_name = "Shepp_reconst_float_64-64-64.raw";
	cond_host.write_detector_name = "detector_float_512-256-180.raw";
	cond_host.write_primary_detector_name = "primary_detector_rectangle_float_512-256-180.raw";
	cond_host.write_efficiency_filter_name = "efficiency_filter_float_512-256.raw";
	cond_host.write_efficiency_correction_name = "efficiency_correction_float_512-256-180.raw";
	cond_host.read_detector_name = "detector_float_512-256-180.raw";
	cond_host.read_primary_detector_name = "primary_detector_rectangle_float_512-256-180.raw";
	cond_host.write_efficiency_map_standardization_name = "efficiency_map_standardization_float_512-256.raw";

	// cond_host.read_fov_name = "fov_5pinhole_5mm_int_512-256.raw";
	// cond_host.read_knife_edge_layer1_name = "layer1_5pinhole_5mm_int_2048-1024.raw";
	// cond_host.read_knife_edge_layer3_name = "layer3_5pinhole_5mm_int_2048-1024.raw";

	cond_host.read_fov_name = "fov_10pinhole_square_2-5mm_int_512-256.raw";
	cond_host.read_knife_edge_layer1_name = "layer1_10pinhole_square_2-5mm_int_2048-1024.raw";
	cond_host.read_knife_edge_layer3_name = "layer3_10pinhole_square_2-5mm_int_2048-1024.raw";

	cond_host.write_origin_fov_name = "origin_fov_int_128-128-20.raw";




	/*---- 時間計測開始&条件表示 start ----*/
  outputLogInit(cond, cond_host);
	/*---- 時間計測開始&条件表示  end -----*/

	// fov作成
	launchCheck(pinhole_theta_xy, pinhole_theta_zx, pinhole_center, cond, cond_host);


	/*----- 時間計測処理 start-----*/
  // outputLogLast(cond, cond_host);
	/*----- 時間計測処理 end-------*/
}

void launchCheck(std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float> pinhole_center, Condition cond, HostCondition cond_host)
{
	// fov作成
	std::vector<int> fov(cond.detector_size_w * cond.detector_size_h, -1);
	std::vector<float> cross_point_collimator_y(cond.pinhole_count * coord_num);
	std::vector<float> cross_point_collimator_z(cond.pinhole_count * coord_num);
	std::vector<float> cross_point_collimator(cond.pinhole_count * coord_num);
	launchCreateGeometry(fov, cross_point_collimator_y, cross_point_collimator_z, pinhole_theta_xy, pinhole_theta_zx, pinhole_center, cond, cond_host);
	launchCheckGeometryCorrect(fov, cross_point_collimator_y, cross_point_collimator_z, cond, cond_host);
}

void launchCreateGeometry(std::vector<int>& fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float> pinhole_center, Condition cond, HostCondition cond_host)
{
	std::vector<int> pinhole_img_layer1(cond.pinhole_img_w * cond.pinhole_img_h, -1);
	std::vector<int> pinhole_img_layer3(pinhole_img_layer1.size(), -1);
	set_pinhole_img(pinhole_img_layer1, pinhole_img_layer3, fov, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, cross_point_collimator_y, cross_point_collimator_z);
	writeRawFile((cond_host.date_directory + cond_host.read_knife_edge_layer1_name).c_str(), pinhole_img_layer1);
	writeRawFile((cond_host.date_directory + cond_host.read_knife_edge_layer3_name).c_str(), pinhole_img_layer3);
	writeRawFile((cond_host.date_directory + cond_host.read_fov_name).c_str(), fov);
}

void launchCheckGeometryCorrect(std::vector<int>& fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, Condition cond, HostCondition cond_host)
{
	std::vector<float> img(cond.img_w * cond.img_h * cond.img_d);
	std::vector<int> origin_fov(cond.img_w * cond.img_h * cond.pinhole_count * 2, -1);
	readRawFile((cond_host.img_file + cond_host.read_activity_map_name).c_str(), img);
	bool is_correct_geometry = CheckGometryCorrect(img, fov, origin_fov, cross_point_collimator_y, cross_point_collimator_z, cond);
	writeRawFile((cond_host.date_directory + cond_host.write_origin_fov_name).c_str(), origin_fov);
	output_test_result(is_correct_geometry);
}

void output_test_result(bool is_correct)
{
	if(is_correct)
	{
		printf("\x1b[32m");     /* 前景色を緑に */
		printf("correct\n");
		printf("\x1b[39m");     /* 前景色をデフォルトに戻す */
	}
	else
	{
		printf("\x1b[31m");     /* 前景色を赤に */
    printf("failed\n");
		printf("\x1b[39m");     /* 前景色をデフォルトに戻す */
	}
}


bool CheckGometryCorrect(std::vector<float>& img, std::vector<int>& fov, std::vector<int>& origin_fov, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z, Condition cond)
{
	rep(pinhole_num, cond.pinhole_count)
	{
		rep(axis, 2)rep(i, cond.detector_size_w)rep(j, cond.detector_size_h)
		{
			if(fov[i * cond.detector_size_w + j] == pinhole_num) continue;
			Eigen::Vector3f cross_point;
			if(axis == 0)
			{
				cross_point(0) = cross_point_collimator_y[coord_num * pinhole_num + X];
				cross_point(1) = cross_point_collimator_y[coord_num * pinhole_num + Y];
				cross_point(2) = cross_point_collimator_y[coord_num * pinhole_num + Z];
			}
			else
			{
				cross_point(0) = cross_point_collimator_z[coord_num * pinhole_num + X];
				cross_point(1) = cross_point_collimator_z[coord_num * pinhole_num + Y];
				cross_point(2) = cross_point_collimator_z[coord_num * pinhole_num + Z];
			}

			Eigen::Vector3f on_detector;
			on_detector(0) = cond.rotation_radius + cond.distance_collimator_to_detector;
			on_detector(1) = transform_origin_coordinate_same_axis(j, cond.detector_size_w, cond.detector_pixel_size_w);
			on_detector(2) = transform_origin_coordinate_opposite_axis(i, cond.detector_size_h, cond.detector_pixel_size_h);

			// cout << "on_detector(X) = " << on_detector(X) << ", on_detector(Y) = " << on_detector(Y) << ", on_detector(Z) = " << on_detector(Z) << endl;
			// cout << "------------------------------------------" << endl;

			Eigen::Vector3f direction_sp = calculate_unit_vector(on_detector, cross_point);
			// 原点上じfovをみる
			showOriginFOV(origin_fov, cross_point, on_detector, direction_sp, cond, pinhole_num, axis);

			//スケーリング処理
			Eigen::Vector3f sp = on_detector / cond.img_pixel_size;

			for (int sample_point = 0; sample_point < 300; sample_point++)
			{
				// cout << "sp(X) = " << sp(X) << ", sp(Y) = " << sp(Y) << ", sp(Z) = " << sp(Z) << endl;

				if(-(cond.img_w - 1.0f) / 2.0f < sp(0) && sp(0) < (cond.img_w - 1.0f) / 2.0f &&
					-(cond.img_h - 1.0f) / 2.0f < sp(1) && sp(1) < (cond.img_h - 1.0f) / 2.0f &&
					-(cond.img_d - 1.0f) / 2.0f < sp(2) && sp(2) < (cond.img_d - 1.0f) / 2.0f )
				{

					int sp_i = transform_img_coordinate_opposite_axis(sp(Y), cond.img_h, 1.);
					int sp_j = transform_img_coordinate_same_axis(sp(X), cond.img_w, 1.);
					int sp_k = transform_img_coordinate_same_axis(sp(Z), cond.img_d, 1.);
					int index = sp_k * cond.img_w * cond.img_h + sp_i * cond.img_w + sp_j;
					if(img[index] > 0.1)
					{
						cout << "------------------------------------------" << endl;
						cout << "img = " << img[index] << endl;
						cout << "pinhole_num = " << pinhole_num << endl;
						cout << "sample_point = " << sample_point << endl;
						cout << "sp(X) = " << sp(X) << ", sp(Y) = " << sp(Y) << ", sp(Z) = " << sp(Z) << endl;
						cout << "i = " << i << ", j = " << j << endl;
						cout << "sp_i = " << sp_i << ", sp_j = " << sp_j << ", sp_k = " << sp_k << endl;
						cout << "axis = " << axis << endl;
						return false;
					}
				}

				sp += direction_sp;
			}
		}
	}
	return true;
}

void showOriginFOV(std::vector<int> &origin_fov, Eigen::Vector3f cross_point, Eigen::Vector3f on_detector, Eigen::Vector3f d, Condition cond, int pinhole_num, int axis)
{
	Eigen::Vector3f on_y_axis;
	on_y_axis(0) = 0.;
	float t = (on_y_axis(0) - on_detector(0)) / d(0);
	on_y_axis(1) = on_detector(1) + t * d(1);
	on_y_axis(2) = on_detector(2) + t * d(2);

	int img_j = transform_img_coordinate_same_axis(on_y_axis(1), cond.img_d, cond.img_pixel_size);
	int img_i = transform_img_coordinate_opposite_axis(on_y_axis(2), cond.img_h, cond.img_pixel_size);

	if(is_out_of_image(img_i, cond.img_h, img_j, cond.img_d)) return;

	origin_fov[(axis + 1) * pinhole_num * cond.img_h * cond.img_d +  cond.img_d * img_i + img_j] = 1;
}


void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, std::vector<int> &fov, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_zx, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z)
{
	for(int pinhole_num = 0; pinhole_num < pinhole_theta_xy.size(); pinhole_num++)
	{
		Eigen::Vector3f center;
		center << pinhole_center[coord_num * pinhole_num + X], pinhole_center[coord_num * pinhole_num + Y], pinhole_center[coord_num * pinhole_num + Z];
		float theta_xy = pinhole_theta_xy[pinhole_num];
		float theta_zx = pinhole_theta_zx[pinhole_num];

		create_pinhole_img_layer1(pinhole_img_layer1, cond, center, theta_xy, theta_zx, pinhole_num, cross_point_collimator_y, cross_point_collimator_z);
		create_pinhole_img_layer3(pinhole_img_layer3, cond, center, theta_xy, theta_zx, pinhole_num, cross_point_collimator_y, cross_point_collimator_z);
		create_fov(fov, cond, center, theta_xy, theta_zx, pinhole_num, cross_point_collimator_y, cross_point_collimator_z);
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

void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 0, cross_point_collimator_y, cross_point_collimator_z); }

void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 1, cross_point_collimator_y, cross_point_collimator_z); }

void create_fov(std::vector<int> &fov, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z)
{ create_pinhole_img( fov, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_zx, pinhole_num, 2, cross_point_collimator_y, cross_point_collimator_z); }


void create_pinhole_img(std::vector<int> &img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_zx, int pinhole_num, int option, std::vector<float>& cross_point_collimator_y, std::vector<float>& cross_point_collimator_z)
{
	/*
	実装のアルゴリズム
	https://paper.dropbox.com/doc/--AjQNF888gr4SxY2nmhuvsAERAQ-zfZOaNRJYJ2QN2vVwrZTU
	*/

	// option { 0: layer1, 1: layer3, 2: fov }
	std::vector<int> width{ cond.pinhole_img_w, cond.pinhole_img_w, cond.detector_size_w };
	std::vector<int> height{ cond.pinhole_img_h, cond.pinhole_img_h, cond.detector_size_h };
	std::vector<float> pixel_size_w{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.detector_pixel_size_w };
	std::vector<float> pixel_size_h{ cond.pinhole_img_pixel_size, cond.pinhole_img_pixel_size, cond.detector_pixel_size_h };
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

	// 交点も持っていく
	if(option == 1)
	{
		cross_point_collimator_y[pinhole_num * coord_num + X] = intersection_edge_y_rot(X);
		cross_point_collimator_y[pinhole_num * coord_num + Y] = intersection_edge_y_rot(Y);
		cross_point_collimator_y[pinhole_num * coord_num + Z] = intersection_edge_y_rot(Z);

		cross_point_collimator_z[pinhole_num * coord_num + X] = intersection_edge_z_rot(X);
		cross_point_collimator_z[pinhole_num * coord_num + Y] = intersection_edge_z_rot(Y);
		cross_point_collimator_z[pinhole_num * coord_num + Z] = intersection_edge_z_rot(Z);
	}


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




std::string makeLogDirectory()
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y_%m_%d_%H.%M.%S/", date);
    std::string s = str;
		s = "log/" + s;

		mkdir("log/", 0777);
		mkdir(s.c_str(), 0777);
		return s;
}


void outputLogLast(Condition cond, HostCondition cond_host)
{
	std::string str;
	std::ostringstream ostr;

 	ostr << "\ntime = " << cond.time << "[s]\n"
			 << "     = " << cond.time / 60 << "[min]\n"
			 << "     = " << cond.time / 3600 << "[h]";

	str = ostr.str();
	std::cout << str << '\n';

	std::string log_text = cond_host.date_directory + "log.txt";
	std::ofstream outputfile(log_text.c_str(), std::ios::app);
	outputfile << str;
	outputfile.close();
}

void outputLogInit(Condition cond, HostCondition cond_host)
{
	std::string str;
	std::ostringstream ostr;
	ostr << "--------------- condition ---------------\n"
			 << "date : " << showCurrentTime() << "\n\n"
			 << "img_w = " << cond.img_w << "\n"
			 << "img_h = " << cond.img_h << "\n"
			 << "img_d = " << cond.img_d << "\n"
			 << "detector_num = " << cond.detector_num << "\n"
			 << "detector_size_w = " << cond.detector_size_w << "\n"
			 << "detector_size_h = " << cond.detector_size_h << "\n"
			 << "surface_source_w = " << cond.surface_source_w << "\n"
			 << "surface_source_h = " << cond.surface_source_h << "\n"
			 << "pinhole_img_w = " << cond.pinhole_img_w << "\n"
			 << "pinhole_img_h = " << cond.pinhole_img_h << "\n"
			 << "pinhole_img_pixel_size = " << cond.pinhole_img_pixel_size << "\n"
			 << "update_count = " << cond.update_count << "\n"
			 << "thread_num = " << cond.thread_num << "\n"
			 << "ray_num = " << cond.ray_num << "\n"
			 << "rotation_radius = " << cond.rotation_radius << "\n"
			 << "distance_collimator_to_detector = " << cond.distance_collimator_to_detector << "\n"
			 << "detector_pixel_size_w = " << cond.detector_pixel_size_w << "\n"
			 << "detector_pixel_size_h = " << cond.detector_pixel_size_h << "\n"
			 << "img_pixel_size = " << cond.img_pixel_size << "\n"
			 << "photon_scale = " << cond.photon_scale << "\n"
			 << "photon_num_efficiency_map = " << cond.photon_num_efficiency_map << "\n"
			 << "max_scattering_count = " << cond.max_scattering_count << "\n"
			 << "cut_off_energy = " << cond.cut_off_energy << "\n"
			 << "medium = " << cond_host.medium << "\n"
			 << "init_position = " << cond_host.init_position << "\n"
			 << "collimator_h = " << cond.collimator_h << "\n";
	if(cond.pinhole_shape == RECTANGLE)
	{
		ostr << "collimator_w = " << cond.collimator_w << "\n"
				 << "aperture_degree = " << cond.aperture_degree << "\n";
	}
	else
	{
		ostr << "collimator_w_y_axis = " << cond.collimator_w_y_axis << "\n"
				 << "collimator_w_z_axis = " << cond.collimator_w_z_axis << "\n"
				 << "aperture_degree_xy = " << cond.aperture_degree_xy << "\n"
				 << "aperture_degree_zx = " << cond.aperture_degree_zx << "\n"
				 << "collimator_w_z_axis = " << cond.collimator_w_z_axis << "\n";
	}

	ostr << "-----------------------------------------\n\n";

	str = ostr.str();
	std::string log_text = cond_host.date_directory + "log.txt";
	std::ofstream outputfile(log_text.c_str(), std::ios::app);
	outputfile << str;
	outputfile.close();

	std::cout << str << std::endl;
}

std::string showCurrentTime()
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y %B %d %A %H:%M:%S", date);
    std::string s = str;
    return s;
}
