/*

doubleにしてある

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)


medium_name
1 : ca
2 : h2o
3 : multi(ca & h2o)

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
// #include "xcom.h"
#define detecter_size 65.

int size = (int)detecter_size;
int photon_num = 100000000;


// -------------------
//medium_num
//1 : ca
//2 : h2o
//3 : multi(ca & h2o)
// -------------------
int medium_num = 1;


float Coherent_ca[201];
float Compton_ca[201];
float Photoelectric_ca[201];
float mu_ca[201];
float Coherent_h2o[201];
float Compton_h2o[201];
float Photoelectric_h2o[201];
float mu_h2o[201];


const char* xcom_csv = "ca_h2o.csv";

void RaySimulation(double* energy_count,double* detecter, float* slab);
void readXcom();


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

	Photon();

	void move();
	void SetProbability(float* slab);
	void ComptonScattering();
};

Photon::Photon() : scatter_(0)
{
	energy_ = 140.;
	float theta = 0.;
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	// ---------------------------------------------
	// 結果の偏りを防ぐため、x座標はランダムに線分で設定
	// -0.0001 <= r <= 0.0001 の範囲内でランダムに設定
	// ---------------------------------------------
	float r = 0.0001 * 2 * (genrand_real1() - 0.5);
	past_(0) = r;
	past_(1) = 0.;
	past_(2) = 0.;

	curr_(0) = r;
	curr_(1) = 0.;
	curr_(2) = 0.;
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

void Photon::SetProbability(float* slab)
{
	float weight = ceil(energy_) - energy_;
	int index1 = floor(energy_);
	int index2 = ceil(energy_);

	// float ca = mu_ca[index1] * weight + mu_ca[index2] * (1. - weight);
	// float h2o = mu_h2o[index1] * weight + mu_h2o[index2] * (1. - weight);
	float ca = mu_ca[index1];
	float h2o = mu_h2o[index1];

	if(medium_num == 1) { mu_max_ = ca; }
	if(medium_num == 2) { mu_max_ = h2o; }
	if(medium_num == 3) { mu_max_ = ca > h2o ? ca : h2o; }

	int j = round(curr_(0) + (128 - 1) / 2.0);
	int i = round((128 - 1) / 2.0 - curr_(1));
	int k = round((128 - 1) - curr_(2));
	if(j >= 128 || i >= 128 || k >= 128 || j < 0 || k < 0 || i < 0)
	{
		printf("j = %d, i = %d, k = %d\n", j, i, k);
		cout << "curr\n" << curr_ << endl;
		cout << "past\n" << past_ << endl;
	}

	medium = slab[k * 128 * 128 + i * 128 + j];

	if(medium > 0.2)
	{
		// coherent_ = Coherent_ca[index1] * weight
		// 			+ Coherent_ca[index2] * (1. - weight);
		//
		// compton_ = Compton_ca[index1] * weight
		// 				+ Compton_ca[index2] * (1. - weight);
		//
		// photo_ = Photoelectric_ca[index1] * weight
		// 				+ Photoelectric_ca[index2] * (1. - weight);

		coherent_ = Coherent_ca[index1];

		compton_ = Compton_ca[index1];

		photo_ = Photoelectric_ca[index1];

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
	// cout << 1. - pow(-1, 2.) << endl;
	init_genrand((unsigned)time(NULL));
	double* energy_count = (double*)calloc(141 * 6, sizeof(double));
	double* detecter = (double*)calloc(size * size * 6, sizeof(double));
	float* slab = (float*)calloc(128 * 128 * 128, sizeof(float));

	string SLAB;

	if(medium_num == 1) 	 { SLAB = "mu-map_ca_float_128-128-128.raw"; }
	else if(medium_num == 2) { SLAB = "mu-map_h2o_float_128-128-128.raw"; }
	else 					 { SLAB = "mu-map_multi_float_128-128-128.raw"; }

	readRawFile(SLAB, 128 * 128 * 128, slab);

	RaySimulation(energy_count, detecter, slab);


	// WriteSeveralFile("h2o", energy_count, detecter);

	string write_csv_name;

	if(medium_num == 1) 	 { write_csv_name = "energy_ca.csv"; }
	else if(medium_num == 2) { write_csv_name = "energy_h2o.csv"; }
	else 					 { write_csv_name = "energy_multi.csv"; }

	for(int i = 0; i < 6; i++)
	{
		double* tmp = (double*)calloc(size * size, sizeof(double));
		for(int j = 0; j < size; j++)
		{
			for(int k = 0; k < size; k++)
				tmp[j * size + k] = detecter[i * size * size + j * size + k];
		}

		char writeFileName[50];
		if(medium_num == 1) 	 { sprintf(writeFileName, "simulation_ca_%02d_double_65-65.raw", i); }
		else if(medium_num == 2) { sprintf(writeFileName, "simulation_h2o_%02d_double_65-65.raw", i); }
		else 					 { sprintf(writeFileName, "simulation_multi_%02d_double_65-65.raw", i); }

		writeRawFile(writeFileName, size * size, tmp);
	}

	FILE* fp;
    if((fp = fopen(write_csv_name.c_str(),"w")) != NULL)
    {
        for(int i = 0; i < 141; i++)
        {
	        fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	i,
	        										energy_count[0 * 141 + i],
					        						energy_count[1 * 141 + i],
					        						energy_count[2 * 141 + i],
					        						energy_count[3 * 141 + i],
					        						energy_count[4 * 141 + i],
					        						energy_count[5 * 141 + i]);
        }
        fclose(fp);
    }
}



void RaySimulation(double* energy_count,double* detecter, float* slab)
{
	int primary_photon_number = 0;
	float* energy_min = (float*)calloc(6, sizeof(float));
	int count_out = 0;
	readXcom();


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
			// cout << i << ": " << mu_ca[i] << "  " << mu_h2o[i] << endl;
		}

		fclose( fp );
	}

	// cout << "mu_ca  mu_h2o" << endl;

	// for(int i = 0; i < 201; i++){ cout << i << ": " << mu_ca << "  " << mu_h2o << endl; }

	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;

    for(int i = 0; i < photon_num; i++)
    {
    	if(i % 1000000 == 0){ printf("photon : %09d\n\n", i); }
    	Photon p;

		p.SetProbability(slab);

		while(1)
		{
			p.move();

			if(abs(p.curr_(0)) > 128 / 2 || abs(p.curr_(1)) > 128 / 2 || p.curr_(2) < 0) { count_out++; break; }

			p.SetProbability(slab);

			if(p.curr_(2) < 10.)
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
					else
					{
						if(p.scatter_ == 0) { break; }
						// p.scatter_++;
					}

					// 散乱回数が5回より大きい＆エネルギーが30より小さい（制限範囲内）
					if(p.scatter_ > 5 || p.energy_ < 30) { break; }
				}
			}
			else // zが10を超えたものの中で検出器を通過したもをのを検出
			{
				float z = 10.;

				if(p.curr_(2) - p.past_(2) < 0.0001) { break; }

				float t = (z - p.past_(2)) / (p.curr_(2) - p.past_(2));
				float x = p.past_(0) + t * (p.curr_(0) - p.past_(0));
				float y = p.past_(1) + t * (p.curr_(1) - p.past_(1));

				if(p.scatter_ == 0) { primary_photon_number++; }

				// xyが検出器の範囲内かどうか
				if(detecter_size / 4 > abs(x) && detecter_size / 4 > abs(y))
				{
					// ピクセルサイズ0.5 cmとするためのスケーリング
					x *= 2;
					y *= 2;

					float j1 = x + (detecter_size) / 2.0;
					float i1 = (detecter_size) / 2.0 - y;

					int J = (int)floor(j1);
					int I = (int)floor(i1);

					int index = (int)round(p.energy_);

					// positionとenergyを検出
					energy_count[p.scatter_ * 141 + index] += 1;
					detecter[p.scatter_ * size * size + I * size + J] += 1;

					if(energy_min[p.scatter_] > p.energy_) { energy_min[p.scatter_] = p.energy_; }
				}
				break;
			}
			p.past_ = p.curr_;
		}
	}

	cout << "primary photon number = " << primary_photon_number << endl;

	if(medium_num == 1)
	{
		cout << "medium : ca" << endl;
		cout << "Difference from the theoretical value : " << 6405584 - primary_photon_number << endl;
	}
	else if(medium_num == 2)
	{
		cout << "medium : h2o" << endl;
		cout << "Difference from the theoretical value : " << 21481029 - primary_photon_number << endl;
	}
	else
		cout << "medium : ca & h2o" << endl;

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
