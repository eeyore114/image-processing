/*

doubleにしてある

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

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
#define DETECTOR_SIZE 65.

const int D_SIZE_INT = (int)DETECTOR_SIZE;
// int photon_num = 100000000;
int photon_num = 100;

// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o) 今回は使わない
//4 : air
// -------------------
int medium_num = 4;


float Coherent_ca[201];
float Compton_ca[201];
float Photoelectric_ca[201];
float mu_ca[201];
float Coherent_h2o[201];
float Compton_h2o[201];
float Photoelectric_h2o[201];
float mu_h2o[201];

const char* XCOM_CSV = "xcomdata_ca_h2o.csv";



void RaySimulation(double* energy_counter,double* detector);
void readXcom();
void WriteEnergySpectrum(double* energy_spectrum, string write_energy_spectrum_name);
void WriteImageProfile(double* detector, string write_image_profile_name);


template <class T>
void writeRawFile (const char fname[], const size_t num, T* image);


class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float energy;
	float theta_cos_, theta_sin_;
	float phi_cos_, phi_sin_;
	float optical_length_;
	int scatter_;
	float coherent_, compton_, photo_, mu_;

	Photon();

	void SetProbability();
	void move();
	void ComptonScattering();
};

Photon::Photon() : scatter_(0)
{
	energy = 140.;
	float theta = 0.;
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);
	past_ << 0., 0., 0.;
}

void Photon::move()
{
	float N_per_N0 = 1. - genrand_real2();
	optical_length_ = - log(N_per_N0) / mu_;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;
}

void Photon::SetProbability()
{
	int index = floor(energy);

	if(medium_num == 1)
	{
		coherent_ = Coherent_ca[index];

		compton_ = Compton_ca[index];

		photo_ = Photoelectric_ca[index];

		mu_ = mu_ca[index];
	}
	else if(medium_num == 2)
	{
		coherent_ = Coherent_h2o[index];

		compton_ = Compton_h2o[index];

		photo_ = Photoelectric_h2o[index];

		mu_ = mu_h2o[index];
	}
	// 空気の場合（move()が動くように適当に設定）
	else { mu_ = 1.; }
}


void Photon::ComptonScattering()
{
	float h = 6.62 * pow(10., -34.);
	float m0 = 9.11 * pow(10., -31.);
	float c = 3.0 * pow(10., 8.);
	float KeV_to_J = 1.602 * pow(10., -19.) * pow(10., 3.);

	energy *= KeV_to_J;

	float lambda = m0 * pow(c, 2.) / energy;
	float relative_frequency = (lambda + 2.) / (9 * lambda + 2);
	float rho;
	int flag = 0;

	while(flag == 0)
	{
		float r1 = genrand_real1();
		float r2 = genrand_real1();
		float r3 = genrand_real1();

		if(relative_frequency > r1)
		{
			rho = 1. + (2. / lambda) * r2;
			flag = (r3 < 4 * (1. / rho - 1 / pow(rho, 2.)));
		}
		else
		{
			rho = (2. + lambda) / (lambda + 2. * (1. - r2));
			flag = (r3 < (pow(lambda - rho * lambda + 1., 2.) + 1 / rho) / 2.);
		}
	}

	float lambda_dash = rho * lambda;
	energy = m0 * pow(c, 2.) / lambda_dash;
	energy /= KeV_to_J;

	float theta_relative_cos = 1. - (lambda_dash - lambda);
	float theta_relative_sin = sqrt(1. - pow(theta_relative_cos, 2.));
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
	double* energy_spectrum = (double*)calloc(141 * 6, sizeof(double));
	double* detector = (double*)calloc(D_SIZE_INT * D_SIZE_INT * 6, sizeof(double));

	if(medium_num == 3) 	 { printf("\x1b[31m"); cout << "multi(ca & h2o) can't select. please change [medium_num]" << endl; printf("\x1b[0m"); exit(-1); }


	RaySimulation(energy_spectrum, detector);

	string write_energy_spectrum_name;
	string write_image_profile_name;

	if(medium_num == 1)
	{
		write_energy_spectrum_name = "result_ca/energy_ca.csv";
		write_image_profile_name = "result_ca/position_ca.csv";
	}
	else if(medium_num == 2)
	{
		write_energy_spectrum_name = "result_h2o/energy_h2o.csv";
		write_image_profile_name = "result_h2o/position_h2o.csv";
	}
	else if (medium_num == 3)
	{
		write_energy_spectrum_name = "result_multi/energy_multi.csv";
		write_image_profile_name = "result_multi/position_multi.csv";
	}
	else
	{
		write_energy_spectrum_name = "result_air/energy_air.csv";
		write_image_profile_name = "result_air/position_air.csv";
	}

	for(int i = 0; i < 6; i++)
	{
		if(medium_num == 4) { break; }

		double* tmp = (double*)calloc(D_SIZE_INT * D_SIZE_INT, sizeof(double));
		for(int j = 0; j < D_SIZE_INT; j++)
		{
			for(int k = 0; k < D_SIZE_INT; k++)
				tmp[j * D_SIZE_INT + k] = detector[i * D_SIZE_INT * D_SIZE_INT + j * D_SIZE_INT + k];
		}

		char writeFileName[50];
		if (medium_num == 1) 	 { sprintf(writeFileName, "result_ca/simulation_ca_%02d_double_65-65.raw", i); }
		else if(medium_num == 2) { sprintf(writeFileName, "result_h2o/simulation_h2o_%02d_double_65-65.raw", i); }

		writeRawFile(writeFileName, D_SIZE_INT * D_SIZE_INT, tmp);
	}

	for(int i = 0; i < 141 * 6; i++) { if(energy_spectrum[i] < 0.01) { energy_spectrum[i] = 1.; } }

	for(int i = 0; i < D_SIZE_INT * D_SIZE_INT * 6; i++) { if(detector[i] < 0.01) { detector[i] = 1.; }	}

	WriteEnergySpectrum(energy_spectrum, write_energy_spectrum_name);
	WriteImageProfile(detector, write_image_profile_name);
}


void RaySimulation(double* energy_counter,double* detector)
{
	int primary_photon_number = 0;
	float* energy_min = (float*)calloc(6, sizeof(float));

	readXcom();

	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;

  for(int i = 0; i < photon_num; i++)
  {
  	if(i % 1000000 == 0){ printf("photon : %09d\n\n", i); }

  	Photon p;

		while(1)
		{
			p.SetProbability();

			p.move();

			if(p.curr_(2) < 10.)
			{
				if(medium_num == 4) { p.past_ = p.curr_; continue; }

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
				if(p.scatter_ > 5 || p.energy < 30) { break; }
			}
			else // zが10を超えたものの中で検出器を通過したもをのを検出
			{
				float z = 10.;

				if(p.curr_(2) - p.past_(2) < 0.0001) { break; }

				float t = (z - p.past_(2)) / (p.curr_(2) - p.past_(2));
				float x = p.past_(0) + t * (p.curr_(0) - p.past_(0));
				float y = p.past_(1) + t * (p.curr_(1) - p.past_(1));

				// xyが検出器の範囲内かどうか
				if(DETECTOR_SIZE / 4 > abs(x) && DETECTOR_SIZE / 4 > abs(y))
				{
					if(p.scatter_ == 0) { primary_photon_number++; }

					// ピクセルサイズ0.5 cmとするためのスケーリング
					x *= 2;
					y *= 2;

					float j1 = x + (DETECTOR_SIZE) / 2.0;
					float i1 = (DETECTOR_SIZE) / 2.0 - y;

					int J = (int)floor(j1);
					int I = (int)floor(i1);

					if(p.scatter_ == 0)
					{
						if(J != 32 || I != 32)
							printf(" I = %d, J = %d\n", I, J);
					}

					int index = (int)round(p.energy);

					// positionとenergyを検出
					energy_counter[p.scatter_ * 141 + index] += 1;
					detector[p.scatter_ * D_SIZE_INT * D_SIZE_INT + I * D_SIZE_INT + J] += 1;

					int test = detector[0 * D_SIZE_INT * D_SIZE_INT + 32 * D_SIZE_INT + 32];

					if(primary_photon_number != test)
					{
						printf("i = %d\n", i);
						printf("x = %f, y = %f, j1 = %f, i1 = %f, J = %d, I = %d, index = %d\n", x, y, j1, i1, J, I, index);
						printf("detector = %d primary_photon_number = %d\n", test, primary_photon_number);
						exit(-1);
					}

					if(energy_min[p.scatter_] > p.energy) { energy_min[p.scatter_] = p.energy; }
				}
				break;
			}
			p.past_ = p.curr_;
		}
	}
	cout << "primary photon number = " << primary_photon_number << endl;

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
			int energy_count = 141;
			fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n", "number", "primary", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5");
			for(int j = 0; j < energy_count; j++)
			{
				fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
													energy_spectrum[0 * energy_count + j],
													energy_spectrum[1 * energy_count + j],
													energy_spectrum[2 * energy_count + j],
													energy_spectrum[3 * energy_count + j],
													energy_spectrum[4 * energy_count + j],
													energy_spectrum[5 * energy_count + j]);
			}
	    fprintf(fp, "\n\n");
    }
    fclose(fp);
}

void WriteImageProfile(double* detector, string write_image_profile_name)
{
	FILE* fp;
	if((fp = fopen(write_image_profile_name.c_str(),"w")) != NULL)
	{
		int position_prof = floor(DETECTOR_SIZE / 2);
		fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n", "position", "primary", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5");
  	for(int j = 0; j < D_SIZE_INT; j++)
  	{
  		fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
    										detector[0 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j],
		        						detector[1 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j],
		        						detector[2 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j],
		        						detector[3 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j],
		        						detector[4 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j],
		        						detector[5 * D_SIZE_INT * D_SIZE_INT + position_prof * D_SIZE_INT + j]);
  	}
    fprintf(fp, "\n\n");
	}
	fclose(fp);
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
