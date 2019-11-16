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
#define DETECTOR_SIZE 180.
#define MUMAP_SIZE 128
#define DETECTOR_NUM 180


const int D_SIZE_INT = (int)DETECTOR_SIZE;


// 1 voxelあたりの光子の数
int photon_num = 1000;
// int photon_num = 100;
// int photon_num = 50;

// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o)
//4 : air
// -------------------
int medium_num = 2;

float Coherent_ca[201];
float Compton_ca[201];
float Photoelectric_ca[201];
float mu_ca[201];
float Coherent_h2o[201];
float Compton_h2o[201];
float Photoelectric_h2o[201];
float mu_h2o[201];

const char* XCOM_CSV = "xcomdata_ca_h2o.csv";

void RaySimulation(double* energy_spectrum,double* detector, float* sphere, float scale_ratio_fantom);
void readXcom();
void WriteEnergySpectrum(double* energy_spectrum, string write_energy_spectrum_name);
void WriteImageProfile(double* detector, string write_image_profile_name);
int JudgeIsair(Eigen::Vector3f curr_, float* sphere, float scale_ratio_fantom);

template <class T>
void readRawFile (string fname, const size_t num, T* image);

template <class T>
void writeRawFile (string fname, const size_t num, T* image);

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

	Photon(int i, int j, int k, float scale_ratio_fantom);

	void move();
	void SetProbability(float* sphere, float scale_ratio_fantom);
	void ComptonScattering();
};

Photon::Photon(int i, int j, int k, float scale_ratio_fantom) : scatter_(0)
{
	energy_ = 140.;
	float rnd = 1. * 2 * (genrand_real1() - 0.5);
	float theta = acos(rnd);
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	float mu_x = - (MUMAP_SIZE - 1.0) / 2.0 + k;
	float mu_y =   (MUMAP_SIZE - 1.0) / 2.0 - j;
	float mu_z =   (MUMAP_SIZE - 1.0) / 2.0 - i;

	mu_x *= scale_ratio_fantom;
	mu_y *= scale_ratio_fantom;
	mu_z *= scale_ratio_fantom;

	// 初期位置にある程度ばらつきをつける
	mu_x += 0.1 * 2 * (genrand_real1() - 0.5);
	mu_y += 0.1 * 2 * (genrand_real1() - 0.5);
	mu_z += 0.1 * 2 * (genrand_real1() - 0.5);

	past_ << mu_x, mu_y, mu_z;
	curr_ << mu_x, mu_y, mu_z;
}

void Photon::move()
{
	float N_per_N0 = 1. - genrand_real2();

	optical_length_ = - log(N_per_N0) / mu_max_;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;
}

void Photon::SetProbability(float* sphere, float scale_ratio_fantom)
{
	int index = floor(energy_);

	float ca = mu_ca[index];
	float h2o = mu_h2o[index];

	if(medium_num == 1) { mu_max_ = ca; }
	if(medium_num == 2) { mu_max_ = h2o; }
	if(medium_num == 3) { mu_max_ = ca > h2o ? ca : h2o; }
	if(medium_num == 4) { mu_max_ = 1; }

	int j = round(curr_(0) / scale_ratio_fantom + (MUMAP_SIZE - 1) / 2.0);
	int i = round((MUMAP_SIZE - 1) / 2.0 - (curr_(1) / scale_ratio_fantom));
	int k = round((MUMAP_SIZE - 1) / 2.0 - (curr_(2) / scale_ratio_fantom));

	if(k > 0 && k < MUMAP_SIZE && i > 0 && i < MUMAP_SIZE && j > 0 && j < MUMAP_SIZE)
	{
		medium = sphere[k * MUMAP_SIZE * MUMAP_SIZE + i * MUMAP_SIZE + j];

		if(medium > 0.2)
		{
			coherent_ = Coherent_ca[index];

			compton_ = Compton_ca[index];

			photo_ = Photoelectric_ca[index];

			mu_ = mu_ca[index];
		}
		else if(medium > 0.1)
		{
			coherent_ = Coherent_h2o[index];

			compton_ = Compton_h2o[index];

			photo_ = Photoelectric_h2o[index];

			mu_ = mu_h2o[index];
		}
		// 空気の場合（move()が動くように適当に設定）
		else { mu_ = 1.; }
	}
}

void Photon::ComptonScattering()
{
	float h = 6.62 * pow(10., -34.);
	float m0 = 9.11 * pow(10., -31.);
	float c = 3.0 * pow(10., 8.);
	float KeV_to_J = 1.602 * pow(10., -19.) * pow(10., 3.);

	energy_ *= KeV_to_J;

	float lambda = m0 * pow(c, 2.) / energy_;
	float relative_frequency = (lambda + 2.) / (9 * lambda + 2);
	float rho;
	int isAcceptable = 0;

	while(isAcceptable == 0)
	{
		float r1 = genrand_real1();
		float r2 = genrand_real1();
		float r3 = genrand_real1();

		if(relative_frequency > r1)
		{
			rho = 1. + (2. / lambda) * r2;
			isAcceptable = (r3 < 4 * (1. / rho - 1 / pow(rho, 2.)));
		}
		else
		{
			rho = (2. + lambda) / (lambda + 2. * (1. - r2));
			isAcceptable = (r3 < (pow(lambda - rho * lambda + 1., 2.) + 1 / rho) / 2.);
		}
	}

	float lambda_dash = rho * lambda;
	energy_ = m0 * pow(c, 2.) / lambda_dash;
	energy_ /= KeV_to_J;

	float theta_relative_cos = 1. - (lambda_dash - lambda);
	float theta_relative_sin = sqrt(1. - pow(theta_relative_cos, 2.));

	if(isnan(theta_relative_sin)) {theta_relative_sin = 0.;}

	// sqrt(0)
	// if(isnan(theta_relative_sin)) { theta_relative_sin = 0.; }

	float phi_relative = genrand_real1() * 2. * M_PI;

	float theta_cos_n = theta_cos_;
	float theta_sin_n = theta_sin_;
	float phi_cos_n = phi_cos_;
	float phi_sin_n = phi_sin_;

	theta_cos_ = -theta_sin_n * theta_relative_sin * cos(phi_relative) + theta_cos_n * theta_relative_cos;
	theta_sin_ = sqrt(1. - pow(theta_cos_, 2.));

	// 0割り対策
	if(theta_sin_ < 0.0000001)
	{
		float rnd_phi = 2 * M_PI * genrand_real2();
		phi_cos_ = cosf(rnd_phi);
		phi_sin_ = sinf(rnd_phi);
	}
	else
	{
		phi_cos_ = (theta_sin_n * phi_cos_n * theta_relative_cos + theta_cos_n * phi_cos_n * theta_relative_sin * cos(phi_relative) - phi_sin_n * theta_relative_sin * sin(phi_relative)) / theta_sin_;
		phi_sin_ = (theta_sin_n * phi_sin_n * theta_relative_cos + theta_cos_n * phi_sin_n * theta_relative_sin * cos(phi_relative) + phi_cos_n * theta_relative_sin * sin(phi_relative)) / theta_sin_;
	}
}

int main()
{
	init_genrand((unsigned)time(NULL));
	double* energy_spectrum = (double*)calloc(4 * 141 * 6, sizeof(double));
	double* detector = (double*)calloc(6 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT, sizeof(double));
	float* sphere = (float*)calloc(MUMAP_SIZE * MUMAP_SIZE * MUMAP_SIZE, sizeof(float));
	// 球ファントムの拡大比
	float scale_ratio_fantom = 0.2;

	string Sphere;

	if(medium_num == 1) 	 { cout << "Ca can't select." << endl; exit(-1); }
	else if(medium_num == 2)
	{
		Sphere = "mu-map_h2o_sphere_float_128-128-128.raw";
		cout << "medium : h2o" << endl;
	}
	else if(medium_num == 3)
	{
		Sphere = "mu-map_multi_sphere_float_128-128-128.raw";
		cout << "medium : h2o & ca" << endl;
	}
	else
	{
		Sphere = "mu-map_air_sphere_float_128-128-128.raw";
		cout << "medium : air" << endl;
	}

	readRawFile(Sphere, MUMAP_SIZE * MUMAP_SIZE * MUMAP_SIZE, sphere);

	RaySimulation(energy_spectrum, detector, sphere, scale_ratio_fantom);

	double* tmp_sum = (double*)calloc(D_SIZE_INT * D_SIZE_INT, sizeof(double));
	for(int i = 0; i < 6; i++)
	{
		double* tmp = (double*)calloc(D_SIZE_INT * D_SIZE_INT, sizeof(double));
		for(int j = 0; j < D_SIZE_INT; j++)
		{
			for(int k = 0; k < D_SIZE_INT; k++)
			{
				tmp[j * D_SIZE_INT + k] += detector[i * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + 0 * D_SIZE_INT * D_SIZE_INT + j * D_SIZE_INT + k];
				tmp_sum[j * D_SIZE_INT + k] += detector[i * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + 0 * D_SIZE_INT * D_SIZE_INT + j * D_SIZE_INT + k];
			}
		}

		char writeFileName[90];
		if(medium_num == 2) { sprintf(writeFileName, "result_h2o/from_voxel/1mm_collimator/sphere_h2o_single_pinhole_%02d_double_180-180.raw", i); }
		else if (medium_num == 3) { sprintf(writeFileName, "result_multi/from_voxel/sphere_multi_single_pinhole_%02d_double_180-180.raw", i); }
		else if (medium_num == 4)  { sprintf(writeFileName, "result_air/from_voxel/sphere_air_single_pinhole_%02d_double_180-180.raw", i); }

		writeRawFile(writeFileName, D_SIZE_INT * D_SIZE_INT, tmp);
	}

	string writeFileName;
	if(medium_num == 2) { writeFileName = "result_h2o/from_voxel/1mm_collimator/sphere_h2o_single_pinhole_sum_double_180-180.raw"; }
	else if (medium_num == 3) { writeFileName = "result_multi/from_voxel/sphere_multi_single_pinhole_sum_double_180-180.raw"; }
	else  { writeFileName = "result_air/from_voxel/sphere_air_single_pinhole_sum_double_180-180.raw"; }
	writeRawFile(writeFileName, D_SIZE_INT * D_SIZE_INT, tmp_sum);

	string write_energy_spectrum_name;
	string write_image_profile_name;

	if(medium_num == 2)
	{
		write_energy_spectrum_name = "result_h2o/from_voxel/1mm_collimator/energy_h2o.csv";
		write_image_profile_name = "result_h2o/from_voxel/1mm_collimator/position_h2o.csv";
	}
	else if (medium_num == 3)
	{
		write_energy_spectrum_name = "result_multi/from_voxel/energy_multi.csv";
		write_image_profile_name = "result_multi/from_voxel/position_multi.csv";
	}
	else
	{
		write_energy_spectrum_name = "result_air/from_voxel/energy_air.csv";
		write_image_profile_name = "result_air/from_voxel/position_air.csv";
	}

	for(int i = 0; i < 4 * 141 * 6; i++) { if(energy_spectrum[i] < 0.01) { energy_spectrum[i] = 1.; } }

	for(int i = 0; i < 6 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT; i++) { if(detector[i] < 0.01) { detector[i] = 1.; }	}

	WriteEnergySpectrum(energy_spectrum, write_energy_spectrum_name);
	WriteImageProfile(detector, write_image_profile_name);

}


void RaySimulation(double* energy_spectrum,double* detector, float* sphere, float scale_ratio_fantom)
{
	// [0] = 0度、[1] = 90度、[2] = 180度、[3] = 270度
	int primary_photon_number[4] = {};
	float* energy_min = (float*)calloc(6, sizeof(float));
	int count_out = 0;
	int count = 0;

	readXcom();

	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;

	for(int i = 0; i < MUMAP_SIZE; i++)
  {
  	printf("processing... %03d / 128\n\n", i);

  	for(int j = 0; j < MUMAP_SIZE; j++)
  	{
  		for(int k = 0; k < MUMAP_SIZE; k++)
  		{
				if(sphere[i * MUMAP_SIZE * MUMAP_SIZE + MUMAP_SIZE * j + k] > 0.0001)
				{
					for(int m = 0; m < photon_num; m++)
					{
			    	Photon p(i, j, k, scale_ratio_fantom);

						p.SetProbability(sphere, scale_ratio_fantom);

						while(1)
						{
							p.move();
							if(isnan(p.curr_.array()).any()) { break; }

							p.SetProbability(sphere, scale_ratio_fantom);

							// if(abs(p.curr_(2)) > 128) { break; }

							int isair = JudgeIsair(p.curr_, sphere, scale_ratio_fantom);

							// if(medium_num == 4) { p.past_ = p.curr_; continue; }

							if(!isair)
							{
								float isDismissal = 1. - genrand_real2();

								if( isDismissal <= p.mu_ / p.mu_max_ )
								{
									float probability = genrand_real2();

									// 光電効果
									if(probability < p.photo_) { break; }
									// コンプトン効果
									else if(probability < p.photo_ + p.compton_)
									{
										p.ComptonScattering();
										p.scatter_++;
									}
									// コヒーレント効果
									else { p.scatter_++; }

									// 散乱回数が5回より大きい＆エネルギーが30より小さい（制限範囲内）
									if(p.scatter_ > 5 || p.energy_ < 30) { break; }
								}
							}
							else
							{
								for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
								{
									float theta = theta_degree * M_PI / 180.;

									// 条件
									float rotation_radius = 10.;
									float distance_collimator_to_detector = 7.5;
									float height_collimator = 1.;
									float width_collimator = 0.1;

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

									float tan_collimator_xy = photon_vec(1) / photon_vec(0);
									float tan_collimator_xz = (curr_rotated(2) - past_rotated(2)) / (curr_rotated(0) - past_rotated(0));

									// コリメータ手前
									Eigen::Vector3f on_colimator;
									on_colimator(0) = rotation_radius - height_collimator / 2.;
									float vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
									on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
									on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
									float colimator_radius = width_collimator / 2. + tan(M_PI / 6);

									if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

									// コリメータ奥
									on_colimator(0) = rotation_radius + height_collimator / 2.;
									vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
									on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
									on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
									colimator_radius = width_collimator / 2. + tan(M_PI / 6);

									if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

									// コリメータ真ん中
									on_colimator(0) = rotation_radius;
									vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
									on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
									on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
									colimator_radius = width_collimator / 2.;

									if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

									//検出器のピクセルサイズが0.5 cmであるため ×2をしている
									// yの値が大きい→検出器の番号は小さい
									float y_on_detector = on_colimator(1) + distance_collimator_to_detector * tan_collimator_xy * 5;
									float z_on_detector = on_colimator(2) + distance_collimator_to_detector * tan_collimator_xz * 5;

									float j0 = D_SIZE_INT / 2. - y_on_detector;
									float i0 = D_SIZE_INT / 2. - z_on_detector;

									if(j0 < 0. || j0 > D_SIZE_INT) { continue; }
									if(i0 < 0. || i0 > D_SIZE_INT) { continue; }

									int j1 = (int)floor(j0);
									int i1 = (int)floor(i0);
									int index = (int)round(p.energy_);

									// positionとenergyを検出
									detector[p.scatter_ * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (theta_degree / 2) * D_SIZE_INT * D_SIZE_INT + i1 * D_SIZE_INT + j1] += 1;

									if(energy_min[p.scatter_] > p.energy_) { energy_min[p.scatter_] = p.energy_; }

									for(int count_result = 0; count_result < 4; count_result++)
									{
										if(theta_degree == count_result * 90)
										{
											energy_spectrum[count_result * 141 * 6 + p.scatter_ * 141 + index] += 1;
											if(p.scatter_ == 0) { primary_photon_number[count_result]++; }
										}
									}
								}
								break;
							}
							p.past_ = p.curr_;
						}
					}
				}
			}
		}
	}

	cout << "----- primary photon number -----" << endl;
	for(int i = 0; i < 4; i++) { printf("%d° = %d\n", i * 90, primary_photon_number[i]); }
	printf("\n");

	if(medium_num == 1) { cout << "medium : ca" << endl; }
	else if(medium_num == 2) { cout << "medium : h2o" << endl; }
	else if(medium_num == 3) { cout << "medium : ca & h2o" << endl; }
	else { cout << "medium : air" << endl; }

	printf("\n");

	for(int i = 0; i < 6; i++)
		cout << "energy_min[" << i << "] = " <<  energy_min[i] << endl;
}

void readXcom()
{
	FILE* fp;
	int ret = 1;
	fp = fopen(XCOM_CSV , "r" );
	if( fp == NULL )
	{
		printf( "failed to open %s\n", XCOM_CSV );
		exit(-1);
	}


	for(int i = 0; i < 201 && ret != EOF; i++)
	{
		fscanf( fp, "%f,%f,%f,%f,%f,%f,%f,%f", &Coherent_ca[i], &Compton_ca[i], &Photoelectric_ca[i], &mu_ca[i], &Coherent_h2o[i], &Compton_h2o[i], &Photoelectric_h2o[i], &mu_h2o[i] );
	}

	fclose( fp );
}


void WriteEnergySpectrum(double* energy_spectrum, string write_energy_spectrum_name)
{
	FILE* fp;
    if((fp = fopen(write_energy_spectrum_name.c_str(),"w")) != NULL)
    {
        for(int i = 0; i < 4; i++)
        {
    		fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n", "number", "primary", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5");
        	for(int j = 0; j < 141; j++)
        	{
        		fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
	        										energy_spectrum[i * 6 * 141 + 0 * 141 + j],
					        						energy_spectrum[i * 6 * 141 + 1 * 141 + j],
					        						energy_spectrum[i * 6 * 141 + 2 * 141 + j],
					        						energy_spectrum[i * 6 * 141 + 3 * 141 + j],
					        						energy_spectrum[i * 6 * 141 + 4 * 141 + j],
					        						energy_spectrum[i * 6 * 141 + 5 * 141 + j]);
        	}
	        fprintf(fp, "\n\n");

        }
        fclose(fp);
    }
}

void WriteImageProfile(double* detector, string write_image_profile_name)
{
	FILE* fp;
    if((fp = fopen(write_image_profile_name.c_str(),"w")) != NULL)
    {
        for(int i = 0; i < 4; i++)
        {
    		fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n", "position", "primary", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5");
        	for(int j = 0; j < D_SIZE_INT; j++)
        	{
        		fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
	        										detector[0 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j],
															detector[1 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j],
															detector[2 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j],
															detector[3 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j],
															detector[4 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j],
															detector[5 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (i * 45) * D_SIZE_INT * D_SIZE_INT + (D_SIZE_INT / 2) * D_SIZE_INT + j]);
        	}
	        fprintf(fp, "\n\n");

        }
        fclose(fp);
    }

}

int JudgeIsair(Eigen::Vector3f curr_, float* sphere, float scale_ratio_fantom)
{
	int j = round(curr_(0) / scale_ratio_fantom + (MUMAP_SIZE - 1) / 2.0);
	int i = round((MUMAP_SIZE - 1) / 2.0 - (curr_(1) / scale_ratio_fantom));
	int k = round((MUMAP_SIZE - 1) / 2.0 - (curr_(2) / scale_ratio_fantom));

	int isair = 0;
	// 空気の場合は1を入れる
	if(k > 0 && k < MUMAP_SIZE && i > 0 && i < MUMAP_SIZE && j > 0 && j < MUMAP_SIZE)
		isair = (0.02 > sphere[k * MUMAP_SIZE * MUMAP_SIZE + i * MUMAP_SIZE + j]);

	return isair;
}


template <class T>
void readRawFile (string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname.c_str());
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}

template <class T>
void writeRawFile(string fname, const size_t num, T* image)
{
	FILE* fp = fopen(fname.c_str(),"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname.c_str());
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname.c_str());
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}
