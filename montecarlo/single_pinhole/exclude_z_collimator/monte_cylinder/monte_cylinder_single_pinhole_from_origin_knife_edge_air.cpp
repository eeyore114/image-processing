/*

doubleにしてある

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)


medium_name
1 : ca
2 : h2o
3 : multi(ca & h2o)


円柱ファントムの拡大比を
float ratio
で設定


*/

#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define detector_size 65.
#define mumap_size 128

int size = (int)detector_size;

// 1 voxelあたりの光子の数
int photon_num = 100000000;
// int photon_num = 10000000;
// int photon_num = 1000000;
// int photon_num = 100000;
// int photon_num = 10000;
// int photon_num = 50;

// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o)
// -------------------
int medium_num = 2;


void RaySimulation(double* energy_spectrum,double* detector, float* cylinder, float scale_ratio_fantom);

class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float energy_;
	float theta_cos_, theta_sin_;
	float phi_cos_, phi_sin_;
	float optical_length_;
	float coherent_, compton_, photo_, mu_, mu_max_;
	int scatter_;
	float medium;
	float phi;

	Photon();

	void move();
};

Photon::Photon() : scatter_(0)
{
	energy_ = 140.;
	float theta = M_PI / 2.;
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	past_ << 0., 0., 0.;
	curr_ << 0., 0., 0.;
}

void Photon::move()
{
	optical_length_ = 1;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;

}

int main()
{
	init_genrand((unsigned)time(NULL));
	double* energy_spectrum = (double*)calloc(4 * 141 * 6, sizeof(double));
	double* detector = (double*)calloc(size * 180 * 6, sizeof(double));
	float* cylinder = (float*)calloc(mumap_size * mumap_size * mumap_size, sizeof(float));
	// 円柱ファントムの拡大比
	float scale_ratio_fantom = 0.2;

	RaySimulation(energy_spectrum, detector, cylinder, scale_ratio_fantom);
}


void RaySimulation(double* energy_spectrum,double* detector, float* cylinder, float scale_ratio_fantom)
{
	// [0] = 0度、[1] = 90度、[2] = 180度、[3] = 270度
	int primary_photon_number[4] = {};
	float* energy_min = (float*)calloc(6, sizeof(float));
	int count_out = 0;
	int count = 0;


	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;


	for(int m = 0; m < photon_num; m++)
	{
    	Photon p;
			if(m % 1000000 == 0) { printf("photon : %d\n", m); }
			if(abs(p.phi) < 0.009999667 || abs(p.phi) > 6.27318564) { count++; }

		while(1)
		{
			p.move();
			if(isnan(p.curr_.array()).any()) { break; }

			if(abs(p.curr_(2)) > 128) { break; }

			if(pow(p.curr_(0), 2.) + pow(p.curr_(1), 2.) > pow(5., 2.))
			{
				for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
				{
					float theta = theta_degree * M_PI / 180.;

					// 条件
					float rotation_radius = 15.;
					float distance_collimator_to_detector = 7.5;
					float height_collimator = 1.;
					float width_collimator = 0.3;

					Eigen::Vector2f past_rotated;
					Eigen::Vector2f curr_rotated;

					past_rotated(0) = p.past_(0) * cos(-theta) - p.past_(1) * sin(-theta);
					past_rotated(1) = p.past_(0) * sin(-theta) + p.past_(1) * cos(-theta);
					curr_rotated(0) = p.curr_(0) * cos(-theta) - p.curr_(1) * sin(-theta);
					curr_rotated(1) = p.curr_(0) * sin(-theta) + p.curr_(1) * cos(-theta);

					if(abs(curr_rotated(0) - past_rotated(0)) < 0.00001) { continue; }

					Eigen::Vector2f photon_vec;
					photon_vec << curr_rotated(0) - past_rotated(0), curr_rotated(1) - past_rotated(1);

					if(abs(photon_vec(0)) < 0.00001) { continue; }
					// if((photon_vec.array() < 0.).all()) { continue; }
					if(photon_vec(0) < 0.) { continue; }

					// float tan_collimator = tan(photon_vec(1) / photon_vec(0));
					float tan_collimator = photon_vec(1) / photon_vec(0);

					// コリメータ手前
					float vec_scale = (rotation_radius - height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					float y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					if(abs(y_on_collimator) > width_collimator / 2. + (height_collimator / 2.) * tan(M_PI / 6)) { continue; }


					// コリメータ真ん中
					vec_scale = (rotation_radius - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					if(abs(y_on_collimator) > width_collimator / 2.) { continue; }

					if(theta_degree == 0 || theta_degree == 359)
					{
							if(abs(p.phi) > 0.009999667 && abs(p.phi) < 6.27318564)
							{
								cout << "p.phi = " << p.phi << "\n";
								cout << "p.past_(0) = " << p.past_(0) << "\n";
								cout << "p.past_(1) = " << p.past_(1) << "\n";
								cout << "p.curr_(0) = " << p.curr_(0) << "\n";
								cout << "p.curr_(1) = " << p.curr_(1) << "\n";
								cout << "past_rotated(0) = " << past_rotated(0) << "\n";
								cout << "past_rotated(1) = " << past_rotated(1) << "\n";
								cout << "curr_rotated(0) = " << curr_rotated(0) << "\n";
								cout << "curr_rotated(1) = " << curr_rotated(1) << "\n";
								cout << "photon_vec(0) = " << photon_vec(0) << "\n";
								cout << "photon_vec(1) = " << photon_vec(1) << "\n";
								cout << "tan_collimator = " << tan_collimator << "\n";
								cout << "y_on_collimator = " << y_on_collimator << "\n";
								exit(-1);
							}
					}


					//検出器のピクセルサイズが0.5 cmであるため ×2をしている
					// yの値が大きい→検出器の番号は小さい
					float y_on_detector = y_on_collimator - distance_collimator_to_detector * tan_collimator * 2;

					float i0 = size / 2. + y_on_detector;

					if(i0 < 0. || i0 > 65.) { continue; }

					int i1 = (int)floor(i0);
					int index = (int)round(p.energy_);

					// positionとenergyを検出
					// detector[p.scatter_ * size * 180 + theta_degree / 2 * size + i1] += 1;

					if(energy_min[p.scatter_] > p.energy_) { energy_min[p.scatter_] = p.energy_; }

					for(int count_result = 0; count_result < 4; count_result++)
					{
						if(theta_degree == count_result * 90)
						{
							// energy_spectrum[count_result * 141 * 6 + p.scatter_ * 141 + index] += 1;
							if(p.scatter_ == 0) { primary_photon_number[count_result]++; }
						}
					}
				}
				break;
			}
			p.past_ = p.curr_;
		}
	}

	cout << "----- primary photon number -----" << endl;
	for(int i = 0; i < 4; i++) { printf("%d° = %d\n", i * 90, primary_photon_number[i]); }
	printf("\n");

	if(medium_num == 1) { cout << "medium : ca" << endl; }
	else if(medium_num == 2) { cout << "medium : h2o" << endl; }
	else { cout << "medium : ca & h2o" << endl; }

	printf("\n");

	printf("count = %d\n", count);
}
