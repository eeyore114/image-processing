/*

floatにしてある

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

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

void RaySimulation();
void readXcom(string read_xcom_name, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o);
void WriteEnergySpectrum(float* energy_spectrum, string write_energy_spectrum_name);
void WriteImageProfile(float* detector, string write_image_profile_name, int detector_size, int detector_num);
int JudgeIsair(Eigen::Vector3f curr_, float* read_img, float scale_ratio_fantom, int mumap_size);
void WriteImage(float* detector, string medium_name, string init_position, int w_collimator_mm, string img_name, int detector_size, int detector_num);



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
	int MUMAP_SIZE;

	Photon(int i, int j, int k, float scale_ratio_fantom, int mumap_size, string init_position);

	void move();
	void SetProbability(float* read_img, string medium_name, float scale_ratio_fantom, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o);
	void ComptonScattering();
};

Photon::Photon(int i, int j, int k, float scale_ratio_fantom, int mumap_size, string init_position) : scatter_(0)
{
	energy_ = 140.;
	float rnd = 1. * 2 * (genrand_real1() - 0.5);
	float theta = acos(rnd);
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);
	MUMAP_SIZE = mumap_size;

	if(init_position == "voxel")
	{
		float mu_x = - (MUMAP_SIZE - 1.0) / 2.0 + k;
		float mu_y =   (MUMAP_SIZE - 1.0) / 2.0 - j;
		float mu_z =   (MUMAP_SIZE - 1.0) / 2.0 - i;

		mu_x /= scale_ratio_fantom;
		mu_y /= scale_ratio_fantom;
		mu_z /= scale_ratio_fantom;

		// 初期位置にある程度ばらつきをつける
		mu_x += 0.1 * 2 * (genrand_real1() - 0.5);
		mu_y += 0.1 * 2 * (genrand_real1() - 0.5);
		mu_z += 0.1 * 2 * (genrand_real1() - 0.5);

		past_ << mu_x, mu_y, mu_z;
		curr_ << mu_x, mu_y, mu_z;
	}
	else if(init_position == "origin")
	{
		past_ << 0., 0., 0.;
		curr_ << 0., 0., 0.;
	}
	else { cout << init_position << " can't select." << endl; exit(-1); }


}

void Photon::move()
{
	float N_per_N0 = 1. - genrand_real2();

	optical_length_ = - log(N_per_N0) / mu_max_;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;
}

void Photon::SetProbability(float* read_img, string medium_name, float scale_ratio_fantom, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o)
{
	int index = floor(energy_);

	float ca = mu_ca[index];
	float h2o = mu_h2o[index];

	if(medium_name == "ca") { mu_max_ = ca; }
	if(medium_name == "h2o") { mu_max_ = h2o; }
	if(medium_name == "multi") { mu_max_ = ca > h2o ? ca : h2o; }
	if(medium_name == "air") { mu_max_ = 1; }

	int j = round(curr_(0) * scale_ratio_fantom + (MUMAP_SIZE - 1) / 2.0);
	int i = round((MUMAP_SIZE - 1) / 2.0 - (curr_(1) * scale_ratio_fantom));
	int k = round((MUMAP_SIZE - 1) / 2.0 - (curr_(2) * scale_ratio_fantom));

	if(k > 0 && k < MUMAP_SIZE && i > 0 && i < MUMAP_SIZE && j > 0 && j < MUMAP_SIZE)
	{
		medium = read_img[k * MUMAP_SIZE * MUMAP_SIZE + i * MUMAP_SIZE + j];

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

int main() { RaySimulation(); }

void RaySimulation()
{
	/*-------------- 条件 --------------*/

	// -----------------------------------------
	// ・基本的に長さの単位は cm
	// ・変数に mm と付いているものは mm
	// -----------------------------------------

	// 1 voxelあたりの光子の数
	// int photon_num = 1E3;
	// int photon_num = 1E10;
	int photon_num = 1E6;

	// -------------------
	//medium_name
	// ca
	// h2o
	// multi(ca & h2o)
	// air
	// -------------------
	string medium_name = "h2o";
	const string read_xcom_name = "xcomdata_ca_h2o.csv";

	float rotation_radius = 5.;
	float distance_collimator_to_detector = 5.;
	float height_collimator = 1.;
	float width_collimator = 0.5;

	// 球ファントムの拡大比
	float scale_ratio_fantom = 5;

	// cylinderとかbrainとか
	string img_name = "sphere";
	// "voxel" or "origin"
	string init_position = "voxel";

	const int MUMAP_SIZE =  32;
	const int DETECTOR_NUM = 180;
	const float DETECTOR_SIZE =  50.;

	/*-------------- 条件ここまで --------------*/

	const int D_SIZE_INT = (int)DETECTOR_SIZE;
	int w_collimator_mm = width_collimator * 10;
	string read_file_name;
	ostringstream ostr;

	// [0] = 0度、[1] = 90度、[2] = 180度、[3] = 270度
	int primary_photon_number[4] = {};
	float* energy_min = (float*)calloc(6, sizeof(float));

	init_genrand((unsigned)time(NULL));
	float* energy_spectrum = (float*)calloc(4 * 141 * 6, sizeof(float));
	float* detector = (float*)calloc(6 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT, sizeof(float));
	float* read_img = (float*)calloc(MUMAP_SIZE * MUMAP_SIZE * MUMAP_SIZE, sizeof(float));

	float Coherent_ca[201];
	float Compton_ca[201];
	float Photoelectric_ca[201];
	float mu_ca[201];
	float Coherent_h2o[201];
	float Compton_h2o[201];
	float Photoelectric_h2o[201];
	float mu_h2o[201];

	readXcom(read_xcom_name, Coherent_ca, Compton_ca, Photoelectric_ca, mu_ca, Coherent_h2o, Compton_h2o, Photoelectric_h2o, mu_h2o);

	for(int i = 0; i < 6; i++) { energy_min[i] = 140.; }

	if(medium_name == "ca") { cout << "Ca can't select." << endl; exit(-1); }
	else
	{
		ostr << "read_img/mu-map_" << medium_name << "_" << img_name << "_float_" << MUMAP_SIZE << "-" << MUMAP_SIZE << "-" << MUMAP_SIZE << ".raw";
		read_file_name = ostr.str();
		ostr.str("");
	}

	readRawFile(read_file_name, MUMAP_SIZE * MUMAP_SIZE * MUMAP_SIZE, read_img);

	for(int i = 0; i < MUMAP_SIZE; i++)
  {
  	if( init_position == "voxel" ) { printf("processing... %03d / %d\n\n", i, MUMAP_SIZE); }

  	for(int j = 0; j < MUMAP_SIZE; j++)
  	{
  		for(int k = 0; k < MUMAP_SIZE; k++)
  		{
				if(read_img[i * MUMAP_SIZE * MUMAP_SIZE + MUMAP_SIZE * j + k] > 0.0001)
				{
					if(init_position == "origin") { i = MUMAP_SIZE; j = MUMAP_SIZE; k = MUMAP_SIZE; }
					for(int m = 0; m < photon_num; m++)
					{
						if(init_position == "origin") { if(m % 1000000 == 0){ printf("photon : %09d\n\n", m); } }
			    	Photon p(i, j, k, scale_ratio_fantom, MUMAP_SIZE, init_position);

						p.SetProbability(read_img, medium_name, scale_ratio_fantom, Coherent_ca, Compton_ca, Photoelectric_ca, mu_ca, Coherent_h2o, Compton_h2o, Photoelectric_h2o, mu_h2o);

						while(1)
						{
							p.move();
							if(isnan(p.curr_.array()).any()) { break; }

							p.SetProbability(read_img, medium_name, scale_ratio_fantom, Coherent_ca, Compton_ca, Photoelectric_ca, mu_ca, Coherent_h2o, Compton_h2o, Photoelectric_h2o, mu_h2o);

							int isair = JudgeIsair(p.curr_, read_img, scale_ratio_fantom, MUMAP_SIZE);

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
								int index_d_num = 360. / DETECTOR_NUM;
								for(int theta_degree = 0; theta_degree < 360; theta_degree += index_d_num)
								{
									float theta = theta_degree * M_PI / 180.;

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
									detector[p.scatter_ * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT + (theta_degree / index_d_num) * D_SIZE_INT * D_SIZE_INT + i1 * D_SIZE_INT + j1] += 1;

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

	cout << "medium : " << medium_name << endl;

	printf("\n");

	for(int i = 0; i < 6; i++) { cout << "energy_min[" << i << "] = " <<  energy_min[i] << endl; }

	WriteImage(detector, medium_name, init_position, w_collimator_mm, img_name, D_SIZE_INT, DETECTOR_NUM);

	// ---------------------
	// csvファイル書き込み
	// ---------------------
	string write_energy_spectrum_name;
	string write_image_profile_name;

	ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << "energy_spectrum_" << medium_name << ".csv";
	write_energy_spectrum_name = ostr.str();
	ostr.str("");

	ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << "profile_" << medium_name << ".csv";
	write_image_profile_name = ostr.str();
	ostr.str("");

	for(int i = 0; i < 4 * 141 * 6; i++) { if(energy_spectrum[i] < 0.01) { energy_spectrum[i] = 1.; } }

	for(int i = 0; i < 6 * DETECTOR_NUM * D_SIZE_INT * D_SIZE_INT; i++) { if(detector[i] < 0.01) { detector[i] = 1.; }	}

	WriteEnergySpectrum(energy_spectrum, write_energy_spectrum_name);
	WriteImageProfile(detector, write_image_profile_name, D_SIZE_INT, DETECTOR_NUM);

}

void readXcom(string read_xcom_name, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o)
{
	FILE* fp;
	int ret = 1;
	fp = fopen(read_xcom_name.c_str() , "r" );
	if( fp == NULL )
	{
		printf( "failed to open %s\n", read_xcom_name.c_str() );
		exit(-1);
	}


	for(int i = 0; i < 201 && ret != EOF; i++)
	{
		fscanf( fp, "%f,%f,%f,%f,%f,%f,%f,%f", &Coherent_ca[i], &Compton_ca[i], &Photoelectric_ca[i], &mu_ca[i], &Coherent_h2o[i], &Compton_h2o[i], &Photoelectric_h2o[i], &mu_h2o[i] );
	}

	fclose( fp );
}


void WriteEnergySpectrum(float* energy_spectrum, string write_energy_spectrum_name)
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

void WriteImageProfile(float* detector, string write_image_profile_name, int detector_size, int detector_num)
{
	FILE* fp;
    if((fp = fopen(write_image_profile_name.c_str(),"w")) != NULL)
    {
        for(int i = 0; i < 4; i++)
        {
    			fprintf(fp,"%s,%s,%s,%s,%s,%s,%s\n", "position", "primary", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5");
        	for(int j = 0; j < detector_size; j++)
        	{
        		fprintf(fp,"%d,%f,%f,%f,%f,%f,%f\n",	j,
	        										detector[0 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j],
															detector[1 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j],
															detector[2 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j],
															detector[3 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j],
															detector[4 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j],
															detector[5 * detector_num * detector_size * detector_size + (i * 45) * detector_size * detector_size + (detector_size / 2) * detector_size + j]);
        	}
	        fprintf(fp, "\n\n");

        }
        fclose(fp);
    }
}

int JudgeIsair(Eigen::Vector3f curr_, float* read_img, float scale_ratio_fantom, int mumap_size)
{
	int j = round(curr_(0) * scale_ratio_fantom + (mumap_size - 1) / 2.0);
	int i = round((mumap_size - 1) / 2.0 - (curr_(1) * scale_ratio_fantom));
	int k = round((mumap_size - 1) / 2.0 - (curr_(2) * scale_ratio_fantom));

	int isair = 1;
	// 空気の場合は1を入れる
	if(k > 0 && k < mumap_size && i > 0 && i < mumap_size && j > 0 && j < mumap_size)
		isair = (0.02 > read_img[k * mumap_size * mumap_size + i * mumap_size + j]);

	return isair;
}


void WriteImage(float* detector, string medium_name, string init_position, int w_collimator_mm, string img_name, int detector_size, int detector_num)
{
	float* tmp_sum = (float*)calloc(detector_size * detector_size, sizeof(float));
	ostringstream ostr;

	for(int i = 0; i < 6; i++)
	{
		float* tmp = (float*)calloc(detector_size * detector_size, sizeof(float));
		for(int j = 0; j < detector_size; j++)
		{
			for(int k = 0; k < detector_size; k++)
			{
				tmp[j * detector_size + k] += detector[i * detector_num * detector_size * detector_size + 0 * detector_size * detector_size + j * detector_size + k];
				tmp_sum[j * detector_size + k] += detector[i * detector_num * detector_size * detector_size + 0 * detector_size * detector_size + j * detector_size + k];
			}
		}

		string write_file_name;

		ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << img_name << "_" << medium_name << "_single_pinhole_" << i << "_float_" <<  detector_size << "-" << detector_size << ".raw";
		write_file_name = ostr.str();
		ostr.str("");
		writeRawFile(write_file_name, detector_size * detector_size, tmp);
	}

	string write_file_name;
	ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << img_name << "_" << medium_name << "_single_pinhole_sum_float_" <<  detector_size << "-" << detector_size << ".raw";
	write_file_name = ostr.str();
	ostr.str("");
	writeRawFile(write_file_name, detector_size * detector_size, tmp_sum);


	float* proj_img = (float*)calloc(detector_num * detector_size * detector_size, sizeof(float));
	float* proj_img_primary = (float*)calloc(detector_num * detector_size * detector_size, sizeof(float));
	for(int i = 0; i < detector_num; i++)
	{
		for(int j = 0; j < detector_size; j++)
		{
			for(int k = 0; k < detector_size; k++)
			{
				proj_img_primary[i * detector_size * detector_size + j * detector_size + k] += detector[0 * detector_num * detector_size * detector_size + i * detector_size * detector_size + j * detector_size + k];
				for(int scatter = 0; scatter < 6; scatter++) { proj_img[i * detector_size * detector_size + j * detector_size + k] += detector[scatter * detector_num * detector_size * detector_size + i * detector_size * detector_size + j * detector_size + k]; }
			}
		}
	}

	/*----- 散乱分全て足し合わせた結果を書き込み -----*/
	string write_proj_img_name;
	ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << "detector_single_pinhole_float_" <<  detector_size << "-" << detector_size << "-" << detector_size << ".raw";
	write_proj_img_name = ostr.str();
	ostr.str("");
	writeRawFile(write_proj_img_name, detector_num * detector_size * detector_size, proj_img);
	/*-------------------------------------------*/

	/*----- primaryの検出結果書き込み -----*/
	string write_proj_img_name_primary;
	ostr << "result_" << medium_name << "/from_" << init_position << "/" << w_collimator_mm << "mm_collimator/" << "primary_single_pinhole_float_" <<  detector_size << "-" << detector_size << "-" << detector_size << ".raw";
	write_proj_img_name_primary = ostr.str();
	ostr.str("");
	writeRawFile(write_proj_img_name_primary, detector_num * detector_size * detector_size, proj_img_primary);
	/*----------------------------------*/
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
