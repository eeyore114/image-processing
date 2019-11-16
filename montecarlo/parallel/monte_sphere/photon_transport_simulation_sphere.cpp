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

// 問題点
// 0, 90度の時にプロファイルがギザギザしちゃう（多分45度も）

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
int photon_num = 1000;
// int photon_num = 50;

// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o)
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

const char* xcom_csv = "xcomdata_ca_h2o.csv";

void RaySimulation(double* energy_spectrum,double* detector, float* sphere, float scale_ratio_fantom);
void readXcom();
void WriteEnergySpectrum(double* energy_spectrum, string write_energy_spectrum_name);
void WriteImageProfile(double* detector, string write_image_profile_name);


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
	float theta = M_PI * genrand_real2();
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	float mu_x = - (mumap_size - 1.0) / 2.0 + k;
	float mu_y =   (mumap_size - 1.0) / 2.0 - j;
	float mu_z =   (mumap_size - 1.0) / 2.0 - i;

	mu_x *= scale_ratio_fantom;
	mu_y *= scale_ratio_fantom;
	mu_z *= scale_ratio_fantom;

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

	if(isnan(curr_.array()).any())
	{
		cout << "curr\n" << curr_ << endl;
		cout << "N_per_N0 = " << N_per_N0 << endl;
		cout << "optical_length_ = " << optical_length_ << endl;
		cout << "theta_sin_ = " << theta_sin_ << endl;
		cout << "theta_cos_ = " << theta_cos_ << endl;
		cout << "phi_sin_ = " << phi_sin_ << endl;
		cout << "phi_cos_ = " << phi_cos_ << endl;
		printf("\n\n");
	}
}

void Photon::SetProbability(float* sphere, float scale_ratio_fantom)
{
	float weight = ceil(energy_) - energy_;
	int index1 = floor(energy_);
	int index2 = ceil(energy_);

	float ca = mu_ca[index1] * weight + mu_ca[index2] * (1. - weight);
	float h2o = mu_h2o[index1] * weight + mu_h2o[index2] * (1. - weight);

	if(medium_num == 1) { mu_max_ = ca; }
	if(medium_num == 2) { mu_max_ = h2o; }
	if(medium_num == 3) { mu_max_ = ca > h2o ? ca : h2o; }

	int j = round(curr_(0) / scale_ratio_fantom + (mumap_size - 1) / 2.0);
	int i = round((mumap_size - 1) / 2.0 - (curr_(1) / scale_ratio_fantom));
	int k = round((mumap_size - 1) / 2.0 - (curr_(2) / scale_ratio_fantom));

	if(k > 0 && k < mumap_size && i > 0 && i < mumap_size && j > 0 && j < mumap_size)
	{
		medium = sphere[k * mumap_size * mumap_size + i * mumap_size + j];

		if(medium > 0.2)
		{
			coherent_ = Coherent_ca[index1] * weight
						+ Coherent_ca[index2] * (1. - weight);

			compton_ = Compton_ca[index1] * weight
							+ Compton_ca[index2] * (1. - weight);

			photo_ = Photoelectric_ca[index1] * weight
							+ Photoelectric_ca[index2] * (1. - weight);

			mu_ = ca;
		}
		else
		{
			coherent_ = Coherent_h2o[index1] * weight
						+ Coherent_h2o[index2] * (1. - weight);

			compton_ = Compton_h2o[index1] * weight
							+ Compton_h2o[index2] * (1. - weight);

			photo_ = Photoelectric_h2o[index1] * weight
							+ Photoelectric_h2o[index2] * (1. - weight);

			mu_ = h2o;
		}
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

	if(isnan(theta_cos_) || isnan(theta_sin_))
	{
		cout << "in the ComptonScattering function" << endl;
		cout << "theta_cos_ = " << theta_cos_ << endl;
		cout << "theta_sin_ = " << theta_sin_ << endl;
		cout << "theta_sin_n = " << theta_sin_n << endl;
		cout << "theta_relative_sin = " << theta_relative_sin << endl;
		cout << "cos(phi_relative) = " << cos(phi_relative) << endl;
		cout << "theta_cos_n = " << theta_cos_n << endl;
		cout << "theta_relative_cos = " << theta_relative_cos << endl;
		printf("\n\n");
	}

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
	double* detector = (double*)calloc(size * 180 * 6, sizeof(double));
	float* sphere = (float*)calloc(mumap_size * mumap_size * mumap_size, sizeof(float));
	// 円柱ファントムの拡大比
	float scale_ratio_fantom = 0.2;


	string Sphere;

	if(medium_num == 1) 	 { cout << "Ca can't select." << endl; exit(-1); }
	else if(medium_num == 2)
	{
		Sphere = "mu-map_h2o_sphere_float_128-128-128.raw";
		cout << "medium : h2o" << endl;
	}
	else
	{
		Sphere = "mu-map_h2o_ca_sphere_float_128-128-128.raw";
		cout << "medium : h2o & ca" << endl;
	}

	readRawFile(Sphere, mumap_size * mumap_size * mumap_size, sphere);

	RaySimulation(energy_spectrum, detector, sphere, scale_ratio_fantom);

	string write_energy_spectrum_name;
	string write_image_profile_name;

	if(medium_num == 1)
	{
		write_energy_spectrum_name = "energy_ca_sphere.csv";
		write_image_profile_name = "position_ca_sphere.csv";
	}
	else if(medium_num == 2)
	{
		write_energy_spectrum_name = "energy_h2o_sphere.csv";
		write_image_profile_name = "position_h2o_sphere.csv";
	}
	else
	{
		write_energy_spectrum_name = "energy_multi_sphere.csv";
		write_image_profile_name = "position_multi_sphere.csv";
	}

	for(int i = 0; i < 4 * 141 * 6; i++) { if(energy_spectrum[i] < 0.01) { energy_spectrum[i] = 1.; } }

	for(int i = 0; i < size * 180 * 6; i++) { if(detector[i] < 0.01) { detector[i] = 1.; }	}

	WriteEnergySpectrum(energy_spectrum, write_energy_spectrum_name);
	WriteImageProfile(detector, write_image_profile_name);

	for(int i = 0; i < 6; i++)
	{
		double* tmp = (double*)calloc(size * size, sizeof(double));
		for(int j = 0; j < 180; j++)
		{
			for(int k = 0; k < size; k++)
				tmp[j * size + k] = detector[i * size * 180 + j * size + k];
		}

		char writeFileName[50];
		if(medium_num == 1) 	 { sprintf(writeFileName, "simulation_sphere_ca_%02d_double_65-180.raw", i); }
		else if(medium_num == 2) { sprintf(writeFileName, "simulation_sphere_h2o_%02d_double_65-180.raw", i); }
		else 					 { sprintf(writeFileName, "simulation_sphere_multi_%02d_double_65-180.raw", i); }

		writeRawFile(writeFileName, size * 180, tmp);

	}


}


void RaySimulation(double* energy_spectrum,double* detector, float* sphere, float scale_ratio_fantom)
{
	int primary_photon_number_0 = 0;
	int primary_photon_number_90 = 0;
	int primary_photon_number_180 = 0;
	int primary_photon_number_270 = 0;
	float* energy_min = (float*)calloc(6, sizeof(float));
	int count_out = 0;
	int count = 0;

	readXcom();

	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;

    for(int i = 0; i < mumap_size; i++)
    {
    	printf("processing... %03d / 128\n\n", i);

    	for(int j = 0; j < mumap_size; j++)
    	{
    		for(int k = 0; k < mumap_size; k++)
    		{
				if(sphere[i * mumap_size * mumap_size + mumap_size * j + k] > 0.1)
				{
					for(int m = 0; m < photon_num; m++)
					{
				    	Photon p(i, j, k, scale_ratio_fantom);

						p.SetProbability(sphere, scale_ratio_fantom);

						while(1)
						{
							p.move();

							p.SetProbability(sphere, scale_ratio_fantom);

							if(abs(p.curr_(2)) > 128) { break; }

							if(sqrt(pow(p.curr_(0), 2.) + pow(p.curr_(1), 2.) + pow(p.curr_(2), 2.)) < 10.)
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
							else // xが10を超えたものの中で検出器を通過したもをのを検出
							{
								Eigen::Vector2f photon_vec;
								photon_vec << p.curr_(0) - p.past_(0), p.curr_(1) - p.past_(1);

								// -------------
								// 投影方向の決定
								// -------------
								float phi_detector;
								int phi_detector_degree_int;
								float r = sqrt(pow(photon_vec(0), 2.) + pow(photon_vec(1), 2.));

								if((photon_vec(0) > 0 && photon_vec(1) > 0) || (photon_vec(0) > 0 && photon_vec(1) < 0)) { phi_detector = atan(photon_vec(1) / photon_vec(0)); }
								else if(photon_vec(0) < 0 && photon_vec(1) > 0) { phi_detector = M_PI - atan(-photon_vec(1) / photon_vec(0)); }
								else if(photon_vec(0) < 0 && photon_vec(1) < 0) { phi_detector = M_PI + atan(photon_vec(1) / photon_vec(0)); }
								else { printf("photon_vec(0) or photon_vec(1) is too small.\n"); }

								float phi_detector_degree = phi_detector * 180. / M_PI;

								if(phi_detector_degree < 0) { phi_detector_degree += 360.;}

								if(int(floor(phi_detector_degree)) % 2 == 0) { phi_detector_degree_int = floor(phi_detector_degree); }
								else										 { phi_detector_degree_int = ceil(phi_detector_degree); }

								if(phi_detector_degree > 359.) { phi_detector_degree_int = 0; }

								float phi_detector_degree_int_to_rad = phi_detector_degree_int * M_PI / 180.;

								// ----------------
								// 円上の座標を求める
								// ----------------

								// y = mx + n とする
								float m = (p.curr_(1) - p.past_(1)) / (p.curr_(0) - p.past_(0));
								float n = p.past_(1) - m * p.past_(0);
								float a = pow(m, 2.) + 1.;
								float b = 2 * m * n;
								float c = pow(n, 2.) - pow(10, 2.);

								// 解の公式
								float cmp_plus = (- b + sqrt(pow(b, 2.) - 4 * a * c)) / (2 * a);
								float cmp_minus = (- b - sqrt(pow(b, 2.) - 4 * a * c)) / (2 * a);

								// warning : ここ解の公式プラスにするかマイナスにするか結構適当に決めてる。
								// float x = cmp_minus;
								float x = abs(p.curr_(0) - cmp_plus) < abs(p.curr_(0) - cmp_minus) ? cmp_plus : cmp_minus;

								float y = m * x + n;

								Eigen::Vector2f photon_vec2;

								photon_vec2(0) = x * cos(-phi_detector_degree_int_to_rad) - y * sin(-phi_detector_degree_int_to_rad);
								photon_vec2(1) = x * sin(-phi_detector_degree_int_to_rad) + y * cos(-phi_detector_degree_int_to_rad);


								// -----------------------------------------
								// xyが検出器の範囲内かどうか
								// スケーリングを考慮して範囲は半分の大きさ
								// -----------------------------------------
								if(detector_size / 4 > abs(photon_vec2(1)))
								{
									// ピクセルサイズ0.5 * 0.5 cm^2とするためのスケーリング
									photon_vec2(1) *= 2;

									float i1 = (detector_size) / 2.0 - photon_vec2(1);

									int I = (int)floor(i1);

									int index = (int)round(p.energy_);

									// positionとenergyを検出
									detector[p.scatter_ * size * 180 + (phi_detector_degree_int / 2) * size + I] += 1;

									if(energy_min[p.scatter_] > p.energy_) { energy_min[p.scatter_] = p.energy_; }

									if(phi_detector_degree_int == 0)
									{
										energy_spectrum[0 * 141 * 6 + p.scatter_ * 141 + index] += 1;
										if(p.scatter_ == 0) { primary_photon_number_0++; }
									}

									if(phi_detector_degree_int == 90)
									{
										energy_spectrum[1 * 141 * 6 + p.scatter_ * 141 + index] += 1;
										if(p.scatter_ == 0) { primary_photon_number_90++; }
									}
									if(phi_detector_degree_int == 180)
									{
										energy_spectrum[2 * 141 * 6 + p.scatter_ * 141 + index] += 1;
										if(p.scatter_ == 0) { primary_photon_number_180++; }
									}
									if(phi_detector_degree_int == 270)
									{
										energy_spectrum[3 * 141 * 6 + p.scatter_ * 141 + index] += 1;
										if(p.scatter_ == 0) { primary_photon_number_270++; }
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
	cout << "0° : " << primary_photon_number_0 << endl;
	cout << "90° : " << primary_photon_number_90 << endl;
	cout << "180° : " << primary_photon_number_180 << endl;
	cout << "270° : " << primary_photon_number_270 << endl;
	printf("\n");

	if(medium_num == 1) { cout << "medium : ca" << endl; }
	else if(medium_num == 2) { cout << "medium : h2o" << endl; }
	else { cout << "medium : ca & h2o" << endl; }

	printf("\n");

	for(int i = 0; i < 6; i++)
		cout << "energy_min[" << i << "] = " <<  energy_min[i] << endl;
}

void readXcom()
{
	FILE* fp;
	int ret = 1;
	fp = fopen(xcom_csv , "r" );
	if( fp == NULL )
	{
		printf( "failed to open %s\n", xcom_csv );
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
        	for(int j = 0; j < size; j++)
        	{
        		fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
	        										detector[0 * 6 * size + i * 90 * size + j],
					        						detector[1 * 6 * size + i * 90 * size + j],
					        						detector[2 * 6 * size + i * 90 * size + j],
					        						detector[3 * 6 * size + i * 90 * size + j],
					        						detector[4 * 6 * size + i * 90 * size + j],
					        						detector[5 * 6 * size + i * 90 * size + j]);
        	}
	        fprintf(fp, "\n\n");

        }
        fclose(fp);
    }

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
