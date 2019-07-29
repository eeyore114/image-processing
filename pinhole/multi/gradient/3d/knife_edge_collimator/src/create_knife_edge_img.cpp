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
	float img_pixel_size;
	float pinhole_img_pixel_size;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float d_width;
	float d_height;
	float time;
	float photon_num;
	float aperture_degree;
} Condition;

void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz);
void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num);
void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num);
void create_pinhole_img(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num, int layer);


int main()
{
	Condition cond;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.pinhole_img_w = 1024 * 2;
	cond.pinhole_img_h = 512 * 2;
	cond.rotation_radius = 13;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.5;
	cond.img_pixel_size = 0.2;
	// cond.pinhole_img_pixel_size = 0.08 / 4;
	cond.pinhole_img_pixel_size = 0.04 / 2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.photon_num = 10000000;
	cond.aperture_degree = 24.;
	std::vector<float> pinhole_theta_xy{ -24., -15., 15., 24., -9., 0., 9., -24., -15., 15., 24. };
	// std::vector<float> pinhole_theta_xy(11, 0.);
	std::vector<float> pinhole_theta_yz{ 	-9.,  -9., -9., -9., 0., 0., 0., 9., 9., 9., 9. };
	// std::vector<float> pinhole_theta_yz(11, 0.);
	float pc_y = - (cond.rotation_radius);
	std::vector<float> pinhole_center
	// { -10., pc_y,  4.2f,
	// 	- 3., pc_y,  4.2f,
	// 	  3., pc_y,  4.2f,
	// 	 10., pc_y,  4.2f,
	// 	-6.2f, pc_y,  0.,
	// 	  0., pc_y,  0.,
	// 	 6.2f, pc_y,  0.,
	// 	-10., pc_y, -4.2f,
 	// 	- 3., pc_y, -4.2f,
 	// 	  3., pc_y, -4.2f,
 	// 	 10., pc_y, -4.2f
	// 									};
	{ -11., pc_y,  4.,
		- 3., pc_y,  4.5f,
		  3., pc_y,  4.5f,
		 11., pc_y,  4.,
		-6.5f, pc_y,  0.,
		  0., pc_y,  0.,
		 6.5f, pc_y,  0.,
		-11., pc_y, -4.,
 		- 3., pc_y, -4.5f,
 		  3., pc_y, -4.5f,
 		 11., pc_y, -4.
										};

	/*---- 時間計測開始&条件表示 start ----*/

	/*---- 時間計測開始&条件表示  end -----*/





	// 変数の定義
	std::vector<int> pinhole_img_layer1(cond.pinhole_img_w * cond.pinhole_img_h, -1);
	std::vector<int> pinhole_img_layer3(pinhole_img_layer1.size(), -1);
	set_pinhole_img(pinhole_img_layer1, pinhole_img_layer3, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_yz);

	// 関数呼び出し
	// launch_test_gradient_pinhole(detector, pinhole_theta_xy, pinhole_theta_yz, pinhole_center, cond);

	/*----- 時間計測処理 start-----*/

	/*----- 時間計測処理 end-------*/
}

void set_pinhole_img(std::vector<int> &pinhole_img_layer1, std::vector<int> &pinhole_img_layer3, Condition cond, std::vector<float> pinhole_center, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz)
{
	for(int pinhole_num = 0; pinhole_num < pinhole_theta_xy.size(); pinhole_num++)
	{
		Eigen::Vector3f center;
		center << pinhole_center[coord_num * pinhole_num + X], pinhole_center[coord_num * pinhole_num + Y], pinhole_center[coord_num * pinhole_num + Z];

		float theta_xy = pinhole_theta_xy[pinhole_num];
		float theta_yz = pinhole_theta_yz[pinhole_num];

		create_pinhole_img_layer1(pinhole_img_layer1, cond, center, theta_xy, theta_yz, pinhole_num);
		create_pinhole_img_layer3(pinhole_img_layer3, cond, center, theta_xy, theta_yz, pinhole_num);
	}
	writeRawFile("../result/layer1_11pinhole_5mm_int_2048-1024.raw", pinhole_img_layer1);
	writeRawFile("../result/layer3_11pinhole_5mm_int_2048-1024.raw", pinhole_img_layer3);
}

void create_pinhole_img_layer1(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_yz, pinhole_num, 1); }

void create_pinhole_img_layer3(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num)
{ create_pinhole_img( pinhole_img, cond, pinhole_center, pinhole_theta_xy, pinhole_theta_yz, pinhole_num, 3); }

void create_pinhole_img(std::vector<int> &pinhole_img, Condition cond, Eigen::Vector3f pinhole_center, float pinhole_theta_xy, float pinhole_theta_yz, int pinhole_num, int layer)
{
	float degree_to_rad = M_PI / 180.;
	float collimator_radius = cond.width_collimator / 2.;
	float aperture = cond.aperture_degree * degree_to_rad;
	Eigen::Vector3f intersection_edge;  intersection_edge << 0, collimator_radius / tan(aperture), 0.;
	if(layer == 1) intersection_edge(1) *= -1.;

	Eigen::Vector3f intersection_edge_rot;
	{
		Eigen::Matrix3f rot_xy, rot_yz;
		Eigen::Vector3f axis; axis << 0, 0, 1;
		rot_xy = Eigen::AngleAxisf( pinhole_theta_xy * degree_to_rad, axis );

		axis << 1, 0, 0;
		rot_yz = Eigen::AngleAxisf( pinhole_theta_yz * degree_to_rad, axis );

		intersection_edge_rot = rot_yz * rot_xy * intersection_edge + pinhole_center;
	}

	float y = (layer == 1) ? (- cond.rotation_radius + cond.height_collimator / 2.) : ( - (cond.rotation_radius + cond.height_collimator / 2.));
	for (int i = 0; i < cond.pinhole_img_h; i++)
	{
		for (int j = 0; j < cond.pinhole_img_w; j++)
		{
			float x = (j - (cond.pinhole_img_w - 1.) / 2.) * cond.pinhole_img_pixel_size;
			// float y = ( - (cond.rotation_radius + cond.height_collimator / 2.));
			float z = ((cond.pinhole_img_h - 1.) / 2. - i) * cond.pinhole_img_pixel_size;
			Eigen::Vector3f img; img << x, y, z;

			Eigen::Vector3f base; base = pinhole_center - intersection_edge_rot;
			Eigen::Vector3f intersection_to_img; intersection_to_img = img - intersection_edge_rot;

			float dot = intersection_to_img.dot(base);
			float norm1 = intersection_to_img.norm();
			float norm2 = base.norm();

			// std::cout << "img" << std::endl;
			// std::cout << img << std::endl;


			if(false)
			{
				std::cout << "------------------" << std::endl;
				std::cout << "pinhole_num : " << pinhole_num << std::endl;
				std::cout << "intersection_edge" << std::endl;
				std::cout << intersection_edge << std::endl;
				std::cout << "intersection_edge_rot" << std::endl;
				std::cout << intersection_edge_rot << std::endl;
				std::cout << "img" << std::endl;
				std::cout << img << std::endl;
				std::cout << "base" << std::endl;
				std::cout << base << std::endl;
				std::cout << "intersection_to_img" << std::endl;
				std::cout << intersection_to_img << std::endl;
				std::cout << "theta = " << acosf( dot / (norm1 * norm2) ) << std::endl;
				std::cout << "aperture = " << aperture << std::endl;
				printf("\n\n\n");
				exit(0);
			}




			if ( acosf( dot / (norm1 * norm2) ) <  aperture ){
				if(false)
				{
					std::cout << "------------------" << std::endl;
					std::cout << "pinhole_num : " << pinhole_num << std::endl;
					std::cout << "intersection_edge" << std::endl;
					std::cout << intersection_edge << std::endl;
					std::cout << "intersection_edge_rot" << std::endl;
					std::cout << intersection_edge_rot << std::endl;
					std::cout << "img" << std::endl;
					std::cout << img << std::endl;
					std::cout << "base" << std::endl;
					std::cout << base << std::endl;
					std::cout << "intersection_to_img" << std::endl;
					std::cout << intersection_to_img << std::endl;
					std::cout << "theta = " << acosf( dot / (norm1 * norm2) ) << std::endl;
					std::cout << "aperture = " << aperture << std::endl;
					printf("\n\n\n");
					// exit(0);
					sleep(1);
				}
				 pinhole_img[cond.pinhole_img_w * i + j] = pinhole_num;  }
			// else printf("not detect\n");
		}
	}
}
