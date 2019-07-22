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
} Condition;

void launch_test_gradient_pinhole(std::vector<float> &detector, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz, std::vector<float> pinhole_center, Condition cond);
void detect_photon(class Photon p, std::vector<float> &detector, Condition cond);
bool passPinhole(class Photon p, float pinhole_theta_xy, float pinhole_theta_yz, Eigen::Vector3f pinhole_center, Condition cond);

class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float theta_;
	float phi_;
	Photon();
	void move();
};

Photon::Photon()
{
	float rnd = 1. * 2 * (genrand_real1() - 0.5f);
	theta_ = acos(rnd);
	phi_ = genrand_real1() * 2 * M_PI;

	past_ << 0., 0., 0.;
	curr_ << 0., 0., 0.;
}

void Photon::move()
{
	float optical_length_ = genrand_real3();

	curr_(0) = past_(0) + optical_length_ * sin(theta_) * cos(phi_);
	curr_(1) = past_(1) + optical_length_ * sin(theta_) * sin(phi_);
	curr_(2) = past_(2) + optical_length_ * cos(theta_);

	if(isnan(curr_.array()).any())
	{
		printf("nan\n\n");
	}
}

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
	launch_test_gradient_pinhole(detector, pinhole_theta_xy, pinhole_theta_yz, pinhole_center, cond);

	/*----- 時間計測処理 start-----*/

	/*----- 時間計測処理 end-------*/
}


void launch_test_gradient_pinhole(std::vector<float> &detector, std::vector<float> pinhole_theta_xy, std::vector<float> pinhole_theta_yz, std::vector<float> pinhole_center, Condition cond)
{
	for(int i = 0; i < cond.photon_num; i++)
	{
		Photon p;
		p.move();

		for(int pinhole_num = 0; pinhole_num < pinhole_theta_xy.size(); pinhole_num++)
		{
			Eigen::Vector3f pin_center;
			int cn = coord_num;
			int pn = pinhole_num;
			pin_center << pinhole_center[cn * pn + X], pinhole_center[cn * pn + Y], pinhole_center[cn * pn + Z];

			bool pass_pinhole = passPinhole(p, pinhole_theta_xy[pinhole_num], pinhole_theta_yz[pinhole_num], pin_center, cond);

			if(!pass_pinhole) continue;
			#ifdef DEBUG
			printf("detect!\n");
			#endif // DEBUG //

			detect_photon(p, detector, cond);
			break;
		}
	}

	writeRawFile("./result/detector_float_512-256.raw", detector);
}

bool passPinhole(class Photon p, float pinhole_theta_xy, float pinhole_theta_yz, Eigen::Vector3f pinhole_center, Condition cond)
{
	// l1の通過判定
	// l3の通過判定
	// l2の通過判定

	#ifdef DEBUG
	std::cout << "pinhole_theta_xy = " << pinhole_theta_xy << '\n';
	std::cout << "pinhole_theta_yz = " << pinhole_theta_yz << '\n';
	#endif // DEBUG //

	Eigen::Vector3f past_rot, curr_rot, axis;
  Eigen::Matrix3f rot_xy, rot_zx;

	// pinhole_center << 0., - (cond.rotation_radius + cond.height_collimator / 2.), 0.;

	/* x軸正方向からy軸正方向に向けての回転を正として回転 */
	axis << 0, 0, 1;
  rot_xy = Eigen::AngleAxisf( - pinhole_theta_xy * M_PI / 180., axis );

	/* z軸正方向からy軸正方向に向けての回転を正として回転 */
  axis << 1, 0, 0;
  rot_zx = Eigen::AngleAxisf( pinhole_theta_yz * M_PI / 180., axis );

	past_rot = rot_zx * rot_xy * (p.past_ - pinhole_center) + pinhole_center;
	curr_rot = rot_zx * rot_xy * (p.curr_ - pinhole_center) + pinhole_center;

	Eigen::Vector3f on_collimator;
	on_collimator(1) = - (cond.rotation_radius + cond.height_collimator / 2.);
	float t = (on_collimator(1) - past_rot(1)) / (curr_rot(1) - past_rot(1));
	on_collimator(0) = past_rot(0) + t * (curr_rot(0) - past_rot(0));
	on_collimator(2) = past_rot(2) + t * (curr_rot(2) - past_rot(2));

	bool is_in_pinhole = pow(on_collimator(0) - pinhole_center(0), 2.) + pow(on_collimator(2) - pinhole_center(2), 2.) < pow(cond.width_collimator / 2., 2.);
	if (!is_in_pinhole)
  {
    #ifdef DEBUG
    // std::cout << "not detect!" << std::endl;
    #endif // DEBUG //
		return false;
  }

	return true;
}

void detect_photon(class Photon p, std::vector<float> &detector, Condition cond)
{
	// 検出処理
	Eigen::Vector3f on_detector;
	on_detector(1) = - (cond.rotation_radius + cond.distance_collimator_to_detector);
	float t = (on_detector(1) - p.past_(1)) / (p.curr_(1) - p.past_(1));
	on_detector(0) = p.past_(0) + t * (p.curr_(0) - p.past_(0));
	on_detector(2) = p.past_(2) + t * (p.curr_(2) - p.past_(2));

	int i = cond.detector_size_h / 2. -  ceilf(on_detector(2) / cond.d_height);
  int j = cond.detector_size_w / 2. + floorf(on_detector(0) / cond.d_width);

	#ifdef DEBUG

	std::cout << " ------------- " << std::endl;
	std::cout << "p.past_ = " << p.past_ << std::endl;
	std::cout << "p.curr_ = " << p.curr_ << std::endl;
	std::cout << "p.theta_ = " << p.theta_ << std::endl;
	std::cout << "p.phi_ = " << p.phi_ << std::endl;
	std::cout << "i = " << i << ", j = " << j << std::endl;
	std::cout << " ------------- " << std::endl;

	#endif // DEBUG //

	if ( 0 > i || i > cond.detector_size_h || 0 > j || j > cond.detector_size_w ) return;

	detector[cond.detector_size_w * i + j] = 1.;

}
