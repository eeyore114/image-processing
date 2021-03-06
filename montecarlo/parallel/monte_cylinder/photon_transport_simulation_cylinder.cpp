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
#define detecter_size 65.

int size = (int)detecter_size;
int photon_num = 100000000;
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


const char* xcom_csv = "ca_h2o.csv";

void RaySimulation(double* energy_count,double* detecter, float* cylinder);
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
	void SetProbability(float* cylinder);
	void ComptonScattering();
};

Photon::Photon() : scatter_(0)
{
	energy_ = 140.;
	float theta = M_PI / 2.;
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	past_ << 0., 0., 0.;
	curr_ << 0., 0., 0.;
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

void Photon::SetProbability(float* cylinder)
{
	float weight = ceil(energy_) - energy_;
	int index1 = floor(energy_);
	int index2 = ceil(energy_);

	float ca = mu_ca[index1] * weight + mu_ca[index2] * (1. - weight);
	float h2o = mu_h2o[index1] * weight + mu_h2o[index2] * (1. - weight);

	if(medium_num == 1) { mu_max_ = ca; }
	if(medium_num == 2) { mu_max_ = h2o; }
	if(medium_num == 3) { mu_max_ = ca > h2o ? ca : h2o; }

	// 円柱ファントムの拡大比
	float ratio = 0.2;

	int j = round(curr_(0) / ratio + (128 - 1) / 2.0);
	int i = round((128 - 1) / 2.0 - (curr_(1) / ratio));
	int k = round((128 - 1) / 2.0 - (curr_(2) / ratio));

	

	//FIXME : 本来なら0のところはkにしなくちゃいけない
	if(k > 0 && k < 128)
	{
		medium = cylinder[0 * 128 * 128 + i * 128 + j];
	}

	

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
	double* detecter = (double*)calloc(size * 180 * 6, sizeof(double));
	float* cylinder = (float*)calloc(128 * 128 * 128, sizeof(float));

	string Cylinder;

	if(medium_num == 1) 	 { cout << "Ca can't select." << endl; exit(-1); } 
	else if(medium_num == 2)
	{
		Cylinder = "mu-map_h2o_cylinder_float_128-128-128.raw";
		cout << "medium : h2o" << endl; 
	} 
	else
	{
		Cylinder = "mu-map_h2o_ca_cylinder_float_128-128-128.raw";
		cout << "medium : h2o & ca" << endl; 
	} 



	readRawFile(Cylinder, 128 * 128 * 128, cylinder);

	RaySimulation(energy_count, detecter, cylinder);



	string write_csv_name;

	if(medium_num == 1) 	 { write_csv_name = "energy_ca_cylinder.csv"; } 
	else if(medium_num == 2) { write_csv_name = "energy_h2o_cylinder.csv"; } 
	else 					 { write_csv_name = "energy_multi_cylinder.csv"; } 


	for(int i = 0; i < 6; i++)
	{
		double* tmp = (double*)calloc(size * size, sizeof(double));
		for(int j = 0; j < 180; j++)
		{
			for(int k = 0; k < size; k++)
				tmp[j * size + k] = detecter[i * size * 180 + j * size + k];
		}

		char writeFileName[50];
		if(medium_num == 1) 	 { sprintf(writeFileName, "simulation_cylinder_ca_%02d_double_65-180.raw", i); } 
		else if(medium_num == 2) { sprintf(writeFileName, "simulation_cylinder_h2o_%02d_double_65-180.raw", i); } 
		else 					 { sprintf(writeFileName, "simulation_cylinder_multi_%02d_double_65-180.raw", i); } 
		
		writeRawFile(writeFileName, size * 180, tmp);
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



void RaySimulation(double* energy_count,double* detecter, float* cylinder)
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

    for(int i = 0; i < photon_num; i++)
    {
    	if(i % 1000000 == 0){ printf("photon : %09d\n\n", i); }
    	Photon p;

		p.SetProbability(cylinder);

		while(1)
		{
			p.move();


			p.SetProbability(cylinder);

			if(abs(p.curr_(2)) > 128) { break; }

			if(sqrt(pow(p.curr_(0), 2.) + pow(p.curr_(1), 2.)) < 10.)
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
				

				// ここの角度設定がうまくいってない
				if((photon_vec(0) > 0 && photon_vec(1) > 0) || (photon_vec(0) > 0 && photon_vec(1) < 0)) { phi_detector = atan(photon_vec(1) / photon_vec(0)); }
				else if(photon_vec(0) < 0 && photon_vec(1) > 0) { phi_detector = M_PI - atan(-photon_vec(1) / photon_vec(0)); }
				else if(photon_vec(0) < 0 && photon_vec(1) < 0) { phi_detector = M_PI + atan(photon_vec(1) / photon_vec(0)); }
				else { printf("photon_vec(0) is too small.\n");	}
				
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
				if(detecter_size / 4 > abs(photon_vec2(1)))
				{
					// ピクセルサイズ0.5 * 0.5 cm^2とするためのスケーリング
					photon_vec2(1) *= 2;

					float i1 = (detecter_size) / 2.0 - photon_vec2(1);

					int I = (int)floor(i1);

					int index = (int)round(p.energy_);
					
					// positionとenergyを検出
					if(phi_detector_degree_int == 0) { energy_count[p.scatter_ * 141 + index] += 1; }
					detecter[p.scatter_ * size * 180 + (phi_detector_degree_int / 2) * size + I] += 1;

					if(energy_min[p.scatter_] > p.energy_) { energy_min[p.scatter_] = p.energy_; }

					if(phi_detector_degree_int == 0) { if(p.scatter_ == 0) { primary_photon_number_0++; } }
					if(phi_detector_degree_int == 90) { if(p.scatter_ == 0) { primary_photon_number_90++; } }
					if(phi_detector_degree_int == 180) { if(p.scatter_ == 0) { primary_photon_number_180++; } }
					if(phi_detector_degree_int == 270) { if(p.scatter_ == 0) { primary_photon_number_270++; } }

				}
				break;
			}
			p.past_ = p.curr_;
		}
	}

	cout << "----- primary photon number -----" << endl;
	cout << "0° : " << primary_photon_number_0 << endl;
	cout << "90° : " << primary_photon_number_90 << endl;
	cout << "180° : " << primary_photon_number_180 << endl;
	cout << "270° : " << primary_photon_number_270 << endl;
	printf("\n");

	if(medium_num == 1)
	{
		cout << "medium : ca" << endl;
	}
	else if(medium_num == 2)
	{
		cout << "medium : h2o" << endl;
		cout << "----- Difference from the theoretical value -----" << endl;
		cout << "0° : " << primary_photon_number_0 - 119339 << endl;
		cout << "90° : " << primary_photon_number_90 - 119339 << endl;
		cout << "180° : " << primary_photon_number_180 - 119339 << endl;
		cout << "270° : " << primary_photon_number_270 - 119339 << endl;
	}
	else
	{
		cout << "medium : ca & h2o" << endl;
		cout << "----- Difference from the theoretical value -----" << endl;
		cout << "0° : " << primary_photon_number_0  - 57712 << endl;
		cout << "90° : " << primary_photon_number_90 - 119339 << endl;
		cout << "180° : " << primary_photon_number_180 - 119339 << endl;
		cout << "270° : " << primary_photon_number_270 - 119339 << endl;
	}

	
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

