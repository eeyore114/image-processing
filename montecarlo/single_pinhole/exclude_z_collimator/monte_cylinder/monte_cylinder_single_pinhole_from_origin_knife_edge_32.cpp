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


5/11 0:07
なんか間違ってる
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
#define detector_size 50.
#define mumap_size 32

int size = (int)detector_size;

// 1 voxelあたりの光子の数
// int photon_num = 100000000;
int photon_num = 10000000;
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

float Coherent_ca[201];
float Compton_ca[201];
float Photoelectric_ca[201];
float mu_ca[201];
float Coherent_h2o[201];
float Compton_h2o[201];
float Photoelectric_h2o[201];
float mu_h2o[201];

const char* xcom_csv = "xcomdata_ca_h2o.csv";

void RaySimulation(double* energy_spectrum,double* detector, float* cylinder, float scale_ratio_fantom);
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

	Photon();

	void move();
	void SetProbability(float* cylinder, float scale_ratio_fantom);
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
		printf("--- in the move function ---\n");
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

void Photon::SetProbability(float* cylinder, float scale_ratio_fantom)
{
	int index = floor(energy_);

	float ca = mu_ca[index];
	float h2o = mu_h2o[index];

	if(medium_num == 1) { mu_max_ = ca; }
	if(medium_num == 2) { mu_max_ = h2o; }
	if(medium_num == 3) { mu_max_ = ca > h2o ? ca : h2o; }

	int j = round(curr_(0) / scale_ratio_fantom + (mumap_size - 1) / 2.0);
	int i = round((mumap_size - 1) / 2.0 - (curr_(1) / scale_ratio_fantom));
	int k = round((mumap_size - 1) / 2.0 - (curr_(2) / scale_ratio_fantom));

	if(k > 0 && k < mumap_size && i > 0 && i < mumap_size && j > 0 && j < mumap_size)
	{
		medium = cylinder[k * mumap_size * mumap_size + i * mumap_size + j];

		if(medium > 0.2)
		{
			coherent_ = Coherent_ca[index];

			compton_ = Compton_ca[index];

			photo_ = Photoelectric_ca[index];

			mu_ = ca;
		}
		else
		{
			coherent_ = Coherent_h2o[index];

			compton_ = Compton_h2o[index];

			photo_ = Photoelectric_h2o[index];

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
	float* cylinder = (float*)calloc(mumap_size * mumap_size * mumap_size, sizeof(float));
	// 円柱ファントムの拡大比
	float scale_ratio_fantom = 0.2;

	string Cylinder;

	if(medium_num == 1) 	 { cout << "Ca can't select." << endl; exit(-1); }
	else if(medium_num == 2)
	{
		Cylinder = "./mkCylinder/mu-map_h2o_cylinder_float_32-32-32.raw";
		cout << "medium : h2o" << endl;
	}
	else
	{
		Cylinder = "mkCylinder/mu-map_h2o_ca_cylinder_float_32-32-32.raw";
		cout << "medium : h2o & ca" << endl;
	}

	readRawFile(Cylinder, mumap_size * mumap_size * mumap_size, cylinder);

	RaySimulation(energy_spectrum, detector, cylinder, scale_ratio_fantom);

	for(int i = 0; i < 6; i++)
	{
		double* tmp = (double*)calloc(size * size, sizeof(double));
		for(int j = 0; j < 180; j++)
		{
			for(int k = 0; k < size; k++)
				tmp[j * size + k] = detector[i * size * 180 + j * size + k];
		}

		char writeFileName[50];
		if(medium_num == 1) 	 { sprintf(writeFileName, "result/cylinder_ca_single_pinhole_from_origin2_%02d_double_50-180.raw", i); }
		else if(medium_num == 2) { sprintf(writeFileName, "result/cylinder_h2o_single_pinhole_from_origin2_%02d_double_50-180.raw", i); }
		else 					 { sprintf(writeFileName, "result/cylinder_multi_single_pinhole_from_origin2_%02d_double_50-180.raw", i); }

		writeRawFile(writeFileName, size * 180, tmp);

	}

	string write_energy_spectrum_name;
	string write_image_profile_name;

	if(medium_num == 1)
	{
		write_energy_spectrum_name = "result/energy_ca_cylinder_single_pinhole_from_origin2.csv";
		write_image_profile_name = "result/position_ca_cylinder_single_pinhole_from_origin2.csv";
	}
	else if(medium_num == 2)
	{
		write_energy_spectrum_name = "result/energy_h2o_cylinder_single_pinhole_from_origin2.csv";
		write_image_profile_name = "result/position_h2o_cylinder_single_pinhole_from_origin2.csv";
	}
	else
	{
		write_energy_spectrum_name = "result/energy_multi_cylinder_single_pinhole_from_origin2.csv";
		write_image_profile_name = "result/position_multi_cylinder_single_pinhole_from_origin2.csv";
	}

	for(int i = 0; i < 4 * 141 * 6; i++) { if(energy_spectrum[i] < 0.01) { energy_spectrum[i] = 1.; } }

	for(int i = 0; i < size * 180 * 6; i++) { if(detector[i] < 0.01) { detector[i] = 1.; }	}

	WriteEnergySpectrum(energy_spectrum, write_energy_spectrum_name);
	WriteImageProfile(detector, write_image_profile_name);


}


void RaySimulation(double* energy_spectrum,double* detector, float* cylinder, float scale_ratio_fantom)
{
	// [0] = 0度、[1] = 90度、[2] = 180度、[3] = 270度
	int primary_photon_number[4] = {};
	float* energy_min = (float*)calloc(6, sizeof(float));
	int count_out = 0;
	int count = 0;

	readXcom();

	for(int i = 0; i < 6; i++)
		energy_min[i] = 140.;


	for(int m = 0; m < photon_num; m++)
	{
    	Photon p;
    	if(m % 1000000 == 0) { printf("photon : %d\n", m); }

		p.SetProbability(cylinder, scale_ratio_fantom);

		while(1)
		{
			p.move();
			if(isnan(p.curr_.array()).any()) { break; }

			p.SetProbability(cylinder, scale_ratio_fantom);

			if(abs(p.curr_(2)) > 128) { break; }

			if(pow(p.curr_(0), 2.) + pow(p.curr_(1), 2.) < pow(2., 2.))
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
			else // 超えた場合
			{
				for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
				{
					float theta = theta_degree * M_PI / 180.;

					// 条件
					float rotation_radius = 5.;
					float distance_collimator_to_detector = 5.;
					float height_collimator = 1.;
					float width_collimator = 0.5;

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
					if(photon_vec(0) < 0.) { continue; }
					float tan_collimator = photon_vec(1) / photon_vec(0);

					// コリメータ手前
					float vec_scale = (rotation_radius - height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					float y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					if(abs(y_on_collimator) > width_collimator / 2. + height_collimator / 2. * tan(M_PI / 6)) { continue; }

					// コリメータ奥
					vec_scale = (rotation_radius + height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					if(abs(y_on_collimator) > width_collimator / 2. + height_collimator / 2. * tan(M_PI / 6)) { continue; }

					// コリメータ真ん中
					vec_scale = (rotation_radius - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
					y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
					if(abs(y_on_collimator) > width_collimator / 2.) { continue; }


					//検出器のピクセルサイズが0.5 cmであるため ×2をしている
					// yの値が大きい→検出器の番号は小さい
					float y_on_detector = y_on_collimator - distance_collimator_to_detector * tan_collimator * 20;

					float i0 = size / 2. + y_on_detector;

					if(i0 < 0. || i0 > detector_size) { continue; }

					int i1 = (int)floor(i0);
					int index = (int)round(p.energy_);

					// positionとenergyを検出
					detector[p.scatter_ * size * 180 + theta_degree / 2 * size + i1] += 1;

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

	cout << "----- primary photon number -----" << endl;
	for(int i = 0; i < 4; i++) { printf("%d° = %d\n", i * 90, primary_photon_number[i]); }
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
	        										detector[0 * 180 * size + i * 45 * size + j],
					        						detector[1 * 180 * size + i * 45 * size + j],
					        						detector[2 * 180 * size + i * 45 * size + j],
					        						detector[3 * 180 * size + i * 45 * size + j],
					        						detector[4 * 180 * size + i * 45 * size + j],
					        						detector[5 * 180 * size + i * 45 * size + j]);
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
