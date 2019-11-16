/*
小さな円のコリメータでやる場合はこれでやればOK
*/

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
#define detector_size 180.
#define mumap_size 128

int size = (int)detector_size;

// 1 voxelあたりの光子の数
// int photon_num = 100000000;
// int photon_num = 10000000;
// int photon_num = 3000000;
int photon_num = 100000;
// int photon_num = 10000;
// int photon_num = 50;

// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o)
// -------------------
int medium_num = 2;



template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);

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
	float rnd = 0.1 * 2 * (genrand_real1() - 0.5);
	float theta = acos(rnd);
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
	double* detector = (double*)calloc(180 * 180, sizeof(double));
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
	double* detector0 = (double*)calloc(180 * 180, sizeof(double));
	double* detector1 = (double*)calloc(180 * 180, sizeof(double));
	double* detector2 = (double*)calloc(180 * 180, sizeof(double));
	double* detector3 = (double*)calloc(180 * 180, sizeof(double));
	double* detector4 = (double*)calloc(180 * 180, sizeof(double));
	double* detector5 = (double*)calloc(180 * 180, sizeof(double));



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
					float width_collimator = 2;

					Eigen::Vector3f past_rotated;
					Eigen::Vector3f curr_rotated;

					past_rotated(0) = p.past_(0) * cos(-theta) - p.past_(1) * sin(-theta);
					past_rotated(1) = p.past_(0) * sin(-theta) + p.past_(1) * cos(-theta);
					past_rotated(2) = p.past_(2);
					curr_rotated(0) = p.curr_(0) * cos(-theta) - p.curr_(1) * sin(-theta);
					curr_rotated(1) = p.curr_(0) * sin(-theta) + p.curr_(1) * cos(-theta);
					curr_rotated(2) = p.curr_(2);

					if(abs(curr_rotated(0) - past_rotated(0)) < 0.00001) { continue; }

					Eigen::Vector2f photon_vec;
					photon_vec << curr_rotated(0) - past_rotated(0), curr_rotated(1) - past_rotated(1);

					if(abs(photon_vec(0)) < 0.00001) { continue; }
					// if((photon_vec.array() < 0.).all()) { continue; }
					if(photon_vec(0) < 0.) { continue; }

					// float tan_collimator = tan(photon_vec(1) / photon_vec(0));
					float tan_collimator_xy = photon_vec(1) / photon_vec(0);
					float tan_collimator_xz = (curr_rotated(2) - past_rotated(2)) / (curr_rotated(0) - past_rotated(0));

					// コリメータ手前
					Eigen::Vector3f on_colimator;
					on_colimator(0) = rotation_radius - height_collimator / 2.;
					float vec_scale = (rotation_radius - height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
					float colimator_radius = width_collimator / 2. + tan(M_PI / 6);


					if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

					// コリメータ真ん中
					on_colimator(0) = rotation_radius;
					vec_scale = (rotation_radius - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
					colimator_radius = width_collimator / 2.;

					if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

					//検出器のピクセルサイズが0.5 cmであるため ×2をしている
					// yの値が大きい→検出器の番号は小さい
					//  FIXME :
					//
					float y_on_detector = on_colimator(1) + distance_collimator_to_detector * tan_collimator_xy * 2;
					float z_on_detector = on_colimator(2) + distance_collimator_to_detector * tan_collimator_xz * 2;

					float j0 = size / 2. - y_on_detector;
					float i0 = size / 2. - z_on_detector;

					if(j0 < 0. || j0 > size) { continue; }
					if(i0 < 0. || i0 > size) { continue; }

					int j1 = (int)floor(j0);
					int i1 = (int)floor(i0);
					int index = (int)round(p.energy_);

					// positionとenergyを検出
					if(p.scatter_ == 0) { detector0[i1 * size + j1] += 1; }
					if(p.scatter_ == 1) { detector1[i1 * size + j1] += 1; }
					if(p.scatter_ == 2) { detector2[i1 * size + j1] += 1; }
					if(p.scatter_ == 3) { detector3[i1 * size + j1] += 1; }
					if(p.scatter_ == 4) { detector4[i1 * size + j1] += 1; }
					if(p.scatter_ == 5) { detector5[i1 * size + j1] += 1; }

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
	writeRawFile("result_air/detector0_double_180-180.raw", size * size, detector0);
	writeRawFile("result_air/detector1_double_180-180.raw", size * size, detector1);
	writeRawFile("result_air/detector2_double_180-180.raw", size * size, detector2);
	writeRawFile("result_air/detector3_double_180-180.raw", size * size, detector3);
	writeRawFile("result_air/detector4_double_180-180.raw", size * size, detector4);
	writeRawFile("result_air/detector5_double_180-180.raw", size * size, detector5);
}


template <class T>
void writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}
