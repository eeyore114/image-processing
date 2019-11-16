/*
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)
*/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <sstream>
#include <sys/stat.h>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#include "fileio.h"

typedef struct {
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int img_w;
	int img_h;
	int img_d;
	int update_count;
	int photon_scale;
	int max_scattering_count;
	float cut_off_energy;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float img_pixel_size;
	float d_width;
	float d_height;
	float time;
	std::string medium;
	std::string init_position;
} Condition;

void RaySimulation(std::vector<float> &activity_map, std::vector<float> &absorp_map, std::vector<float> &detector, Condition cond);
void readXcom(std::string read_xcom_name, std::vector<float> &Coherent_ca, std::vector<float> &Compton_ca, std::vector<float> &Photoelectric_ca, std::vector<float> &mu_ca,	std::vector<float> &Coherent_h2o, std::vector<float> &Compton_h2o, std::vector<float> &Photoelectric_h2o, std::vector<float> &mu_h2o);
bool JudgeIsair(Eigen::Vector3f curr_, std::vector<float> &absorp_map, Condition cond);
void WriteImage(std::vector<float> &detector, Condition cond);
void outputLog(Condition cond, bool is_console = true, bool is_finish = false);
void make_directory(Condition cond, float w_collimator_mm);
std::string showCurrentTime();
std::string makeLogDirectory();


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

	Photon(int i, int j, int k, Condition cond);

	void move();
	void SetProbability(std::vector<float> &absorp_map, Condition cond, std::vector<float> &Coherent_ca, std::vector<float> &Compton_ca, std::vector<float> &Photoelectric_ca, std::vector<float> &mu_ca,	std::vector<float> &Coherent_h2o, std::vector<float> &Compton_h2o, std::vector<float> &Photoelectric_h2o, std::vector<float> &mu_h2o);
	void ComptonScattering();
};

Photon::Photon(int i, int j, int k, Condition cond) : scatter_(0)
{
	energy_ = 140.;
	float rnd = 1. * 2 * (genrand_real1() - 0.5);
	float theta = acos(rnd);
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);

	if(cond.init_position == "voxel")
	{
		float mu_x = - (cond.img_w - 1.0) / 2.0 + k;
		float mu_y =   (cond.img_h - 1.0) / 2.0 - j;
		float mu_z =   (cond.img_d - 1.0) / 2.0 - i;

		mu_x *= cond.img_pixel_size;
		mu_y *= cond.img_pixel_size;
		mu_z *= cond.img_pixel_size;

		// 初期位置にある程度ばらつきをつける
		mu_x += 0.1 * 2 * (genrand_real1() - 0.5);
		mu_y += 0.1 * 2 * (genrand_real1() - 0.5);
		mu_z += 0.1 * 2 * (genrand_real1() - 0.5);

		past_ << mu_x, mu_y, mu_z;
		curr_ << mu_x, mu_y, mu_z;
	}
	else if(cond.init_position == "origin")
	{
		past_ << 0., 0., 0.;
		curr_ << 0., 0., 0.;
	}
	else { std::cout << cond.init_position << " can't select." << std::endl; exit(-1); }
}

void Photon::move()
{
	float N_per_N0 = 1. - genrand_real2();

	optical_length_ = - log(N_per_N0) / mu_max_;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;
}

void Photon::SetProbability(std::vector<float> &absorp_map, Condition cond, std::vector<float> &Coherent_ca, std::vector<float> &Compton_ca, std::vector<float> &Photoelectric_ca, std::vector<float> &mu_ca,	std::vector<float> &Coherent_h2o, std::vector<float> &Compton_h2o, std::vector<float> &Photoelectric_h2o, std::vector<float> &mu_h2o)
{
	int index = floor(energy_);

	float ca = mu_ca[index];
	float h2o = mu_h2o[index];
 	mu_max_ = ca > h2o ? ca : h2o;

	int j = round(curr_(0) / cond.img_pixel_size + (cond.img_w - 1) / 2.0);
	int i = round((cond.img_h - 1) / 2.0 - (curr_(1) / cond.img_pixel_size));
	int k = round((cond.img_d - 1) / 2.0 - (curr_(2) / cond.img_pixel_size));

	if(k > 0 && k < cond.img_d && i > 0 && i < cond.img_h && j > 0 && j < cond.img_w)
	{
		medium = absorp_map[k * cond.img_w + cond.img_h + i * cond.img_w + j];

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

int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.img_d = 128;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 10;
	cond.photon_scale = 1;
	cond.max_scattering_count = 5;
	cond.cut_off_energy = 30.;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.2;
	cond.update_count = 30;
	cond.img_pixel_size = 0.1;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.medium = "Shepp"; //"Shepp", "Brain", "Sphere"　など（大文字始まり）
	cond.init_position = "voxel"; //"voxel" or "origin"

	/*---- 時間計測開始&条件表示 start ----*/
	clock_t start = clock();
	outputLog(cond);
	/*---- 時間計測開始&条件表示  end -----*/

	std::vector<float> activity_map(cond.img_w * cond.img_h * cond.img_d, 0.);
	std::vector<float> absorp_map(activity_map.size(), 0.);
	int scatter_num = 6;
	std::vector<float> detector(scatter_num * cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);

	std::string read_activity_map_name;
	std::ostringstream ostr;
	ostr << "read_img/" << cond.medium << "_float_" << cond.img_w << "-" << cond.img_h << "-" << cond.img_h << ".raw";
	read_activity_map_name = ostr.str();
	ostr.str("");
	readRawFile(read_activity_map_name, activity_map);

	std::string read_absorp_map_name;
	ostr << "read_img/" << cond.medium << "_absorp_map" << "_float_" << cond.img_w << "-" << cond.img_h << "-" << cond.img_h << ".raw";
	read_absorp_map_name = ostr.str();
	ostr.str("");
	readRawFile(read_absorp_map_name, absorp_map);

	RaySimulation(activity_map, absorp_map, detector, cond);

	int w_collimator_mm = cond.width_collimator * 10;
	/*---- 時間計測終了&ログ書き込み start ----*/
	clock_t end = clock();
  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
	bool is_console = false;
	bool is_finish = true;
  outputLog(cond, is_console, is_finish);

  std::cout << "time = " << cond.time << "[s]\n"
  					 <<	"     = " << cond.time / 60 << "[min]\n"
  					 << "     = " << cond.time / 3600 << "[h]" << std::endl;
	/*---- 時間計測終了&ログ書き込み  end  ----*/


}

void RaySimulation(std::vector<float> &activity_map, std::vector<float> &absorp_map, std::vector<float> &detector, Condition cond)
{
	init_genrand((unsigned)time(NULL));

	std::vector<float> mu_ca(201, 0.);
	std::vector<float> coherent_ca(mu_ca.size(), 0.);
	std::vector<float> compton_ca(mu_ca.size(), 0.);
	std::vector<float> photoelectric_ca(mu_ca.size(), 0.);
	std::vector<float> mu_h2o(mu_ca.size(), 0.);
	std::vector<float> coherent_h2o(mu_ca.size(), 0.);
	std::vector<float> compton_h2o(mu_ca.size(), 0.);
	std::vector<float> photoelectric_h2o(mu_ca.size(), 0.);

	std::string read_xcom_name = "xcomdata_ca_h2o.csv";
	readXcom(read_xcom_name, coherent_ca, compton_ca, photoelectric_ca, mu_ca, coherent_h2o, compton_h2o, photoelectric_h2o, mu_h2o);

	for(int i = 0; i < cond.img_h; i++)
  {
  	if( cond.init_position == "voxel" ) { printf("processing... %03d / %d\n\n", i, cond.img_h); }

  	for(int j = 0; j < cond.img_w; j++)
  	{
  		for(int k = 0; k < cond.img_d; k++)
  		{
				int array_count = k * cond.img_w * cond.img_h + cond.img_w * i + j;

				for(int m = 0; m < activity_map[array_count] * cond.photon_scale; m++)
				{
		    	Photon p(i, j, k, cond);

					p.SetProbability(absorp_map, cond, coherent_ca, compton_ca, photoelectric_ca, mu_ca, coherent_h2o, compton_h2o, photoelectric_h2o, mu_h2o);

					while(true)
					{
						p.move();
						if(isnan(p.curr_.array()).any()) { break; }

						p.SetProbability(absorp_map, cond, coherent_ca, compton_ca, photoelectric_ca, mu_ca, coherent_h2o, compton_h2o, photoelectric_h2o, mu_h2o);

						bool isair = JudgeIsair(p.curr_, absorp_map, cond);

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
								if(p.scatter_ > cond.max_scattering_count || p.energy_ < cond.cut_off_energy) { break; }
							}
						}
						else
						{
							int delta_detector = 360. / cond.detector_num;
							for(int theta_degree = 0; theta_degree < 360; theta_degree += delta_detector)
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
								on_colimator(0) = cond.rotation_radius - cond.height_collimator / 2.;
								float vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
								on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
								on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
								float colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

								if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

								// コリメータ奥
								on_colimator(0) = cond.rotation_radius + cond.height_collimator / 2.;
								vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
								on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
								on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
								colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

								if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

								// コリメータ真ん中
								on_colimator(0) = cond.rotation_radius;
								vec_scale = (on_colimator(0) - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
								on_colimator(1) = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
								on_colimator(2) = past_rotated(2) + vec_scale * (curr_rotated(2) - past_rotated(2));
								colimator_radius = cond.width_collimator / 2.;

								if(pow(on_colimator(1), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { continue; }

								//検出器のピクセルサイズが0.5 cmであるため ×2をしている
								// yの値が大きい→検出器の番号は小さい
								float y_on_detector = on_colimator(1) + cond.distance_collimator_to_detector * tan_collimator_xy * 5;
								float z_on_detector = on_colimator(2) + cond.distance_collimator_to_detector * tan_collimator_xz * 5;

								float j0 = cond.detector_size_w / 2. - y_on_detector;
								float i0 = cond.detector_size_h / 2. - z_on_detector;

								if(j0 < 0. || j0 > cond.detector_size_w) { continue; }
								if(i0 < 0. || i0 > cond.detector_size_h) { continue; }

								int j1 = (int)floor(j0);
								int i1 = (int)floor(i0);
								int index = (int)round(p.energy_);

								// positionとenergyを検出
								detector[p.scatter_ * cond.detector_num * cond.detector_size_w * cond.detector_size_h + (theta_degree / delta_detector) * cond.detector_size_w * cond.detector_size_h + i1 * cond.detector_size_w + j1] += 1;
							}
							break;
						}
						p.past_ = p.curr_;
					}
				}
			}
		}
	}
	WriteImage(detector, cond);
}

void readXcom(std::string read_xcom_name, std::vector<float> &Coherent_ca, std::vector<float> &Compton_ca, std::vector<float> &Photoelectric_ca, std::vector<float> &mu_ca,	std::vector<float> &Coherent_h2o, std::vector<float> &Compton_h2o, std::vector<float> &Photoelectric_h2o, std::vector<float> &mu_h2o)
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
bool JudgeIsair(Eigen::Vector3f curr_, std::vector<float> &absorp_map, Condition cond)
{
	int j = round(curr_(0) / cond.img_pixel_size + (cond.img_w - 1) / 2.0);
	int i = round((cond.img_h - 1) / 2.0 - (curr_(1) / cond.img_pixel_size));
	int k = round((cond.img_d - 1) / 2.0 - (curr_(2) / cond.img_pixel_size));

	bool isair = true;
	// 空気の場合は1を入れる
	if(k > 0 && k < cond.img_d && i > 0 && i < cond.img_h && j > 0 && j < cond.img_w)
		isair = (0.02 > absorp_map[k * cond.img_w * cond.img_h + i * cond.img_w + j]);

	return isair;
}

void WriteImage(std::vector<float> &detector, Condition cond)
{
	std::ostringstream ostr;

	std::vector<float> proj_img_sum(cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> proj_img_primary(proj_img_sum.size(), 0.);
	for(int i = 0; i < cond.detector_num; i++)
	{
		for(int j = 0; j < cond.detector_size_h; j++)
		{
			for(int k = 0; k < cond.detector_size_w; k++)
			{
				proj_img_primary[i * cond.detector_size_w * cond.detector_size_h + j * cond.detector_size_w + k] += detector[0 * cond.detector_num * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w * cond.detector_size_h + j * cond.detector_size_w + k];
				for(int scatter = 0; scatter <= cond.max_scattering_count ; scatter++) { proj_img_sum[i * cond.detector_size_w * cond.detector_size_h + j * cond.detector_size_w + k] += detector[scatter * cond.detector_num * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w * cond.detector_size_h + j * cond.detector_size_w + k]; }
			}
		}
	}

	float w_collimator_mm = cond.width_collimator * 10;
	make_directory(cond, w_collimator_mm);
	/*----- 散乱分全て足し合わせた結果を書き込み -----*/
	std::string write_proj_img_name;
	ostr << "result_" << cond.medium << "/from_" << cond.init_position << "/" << w_collimator_mm << "mm_collimator/" << "detector_single_pinhole_float_" << cond.detector_size_w << "-" << cond.detector_size_h << "-" << cond.detector_num << ".raw";
	write_proj_img_name = ostr.str();
	ostr.str("");
	writeRawFile(write_proj_img_name, proj_img_sum);
	/*-------------------------------------------*/

	/*----- primaryの検出結果書き込み -----*/
	std::string write_proj_img_name_primary;
	ostr << "result_" << cond.medium << "/from_" << cond.init_position << "/" << w_collimator_mm << "mm_collimator/" << "primary_single_pinhole_float_" << cond.detector_size_w << "-" << cond.detector_size_h << "-" << cond.detector_num << ".raw";
	write_proj_img_name_primary = ostr.str();
	ostr.str("");
	writeRawFile(write_proj_img_name_primary, proj_img_primary);
	/*----------------------------------*/
}

void make_directory(Condition cond, float w_collimator_mm)
{
	std::string path;
	std::ostringstream ostr;
	ostr << "result_" << cond.medium << "/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
	ostr << "from_" << cond.init_position << "/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
	ostr << w_collimator_mm << "mm_collimator/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
}

void outputLog(Condition cond, bool is_console, bool is_finish)
{
	std::string str;
	std::ostringstream ostr;
	ostr << "--------------- condition ---------------\n"
			 << "date : " << showCurrentTime() << "\n\n"
			 << "img_w = " << cond.img_w << "\n"
			 << "img_h = " << cond.img_h << "\n"
			 << "img_d = " << cond.img_d << "\n"
			 << "detector_num = " << cond.detector_num << "\n"
			 << "detector_size_w = " << cond.detector_size_w << "\n"
			 << "detector_size_h = " << cond.detector_size_h << "\n"
			 << "update_count = " << cond.update_count << "\n"
			 << "rotation_radius = " << cond.rotation_radius << "\n"
			 << "distance_collimator_to_detector = " << cond.distance_collimator_to_detector << "\n"
			 << "height_collimator = " << cond.height_collimator << "\n"
			 << "width_collimator = " << cond.width_collimator << "\n"
			 << "d_width = " << cond.d_width << "\n"
			 << "d_height = " << cond.d_height << "\n"
			 << "img_pixel_size = " << cond.img_pixel_size << "\n"
			 << "photon_scale = " << cond.photon_scale << "\n"
			 << "max_scattering_count = " << cond.max_scattering_count << "\n"
			 << "cut_off_energy = " << cond.cut_off_energy << "\n"
			 << "medium = " << cond.medium << "\n"
			 << "init_position = " << cond.init_position << "\n"
			 << "-----------------------------------------\n\n";

	if(is_finish)
 	{
 		ostr << "\ntime = " << cond.time << "[s]\n"
			 	 << "     = " << cond.time / 60 << "[min]\n"
			   << "     = " << cond.time / 3600 << "[h]";
 	}

	str = ostr.str();

	if(is_console) { std::cout << str << std::endl;	}
	else
	{
		std::string log_text = makeLogDirectory() + "log.txt";
		std::ofstream outputfile(log_text);
		outputfile << str;
		outputfile.close();
	}
}

std::string showCurrentTime()
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y %B %d %A %p %I:%M:%S", date);
    std::string s = str;
    return s;
}

std::string makeLogDirectory()
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y_%m_%d_%H.%M.%S/", date);
    std::string s = str;
		s = "log/" + s;

		mkdir("log/", 0777);
		mkdir(s.c_str(), 0777);
		return s;
}
