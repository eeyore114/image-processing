/*
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

・多分大丈夫

・条件揃える
・cudaで実装
これを使ってMLEM



*/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"


void SensitivityCorrection(float* f, float* g, struct Condition cond);
void MakeSensitivityMap(float* f, struct Condition cond);

struct Condition {
	int photon_num;
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int surface_source_w;
	int surface_source_h;
	float distance_collimator_to_detector;
	float pixel_size_detector;
	float collimator_h;
	float collimator_w;
};

class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float theta;
	float phi;

	Photon(int i, int j, struct Condition cond);

};

Photon::Photon(int i, int j, struct Condition cond)
{
	float rnd = 1. * 2 * (genrand_real1() - 0.5);
	theta = acos(rnd);
	// theta = 0.;
	// phi = genrand_real1() * 2 * M_PI;
	phi = (genrand_real1() - 0.5) * M_PI;

	float mu_x = - (cond.surface_source_w - 1.0) / 2.0 + j;
	float mu_y =   0.;
	float mu_z =   (cond.surface_source_h - 1.0) / 2.0 - i;

	// float mu_x = 0.;
	// float mu_y = 0.;
	// float mu_z = 0.;

	// 初期位置にある程度ばらつきをつける
	// mu_x += 0.1 * 2 * (genrand_real1() - 0.5);
	// mu_z += 0.1 * 2 * (genrand_real1() - 0.5);

	mu_x += (genrand_real1() - 0.5);
	mu_z += (genrand_real1() - 0.5);

	past_ << mu_x, mu_y, mu_z;
	curr_ << mu_x, mu_y, mu_z;
}



template <class T>
void readRawFile (std::string fname, const size_t num, T* image);

template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);

void mkEfficiencyMap(float* surface_source, float* efficiency_map, struct Condition cond);

int canPassCollimator(class Photon p, Eigen::Vector3f on_colimator, int theta_degree, struct Condition cond);



int main()
{
	/*
	・光子を出す位置を判定する画像を読み込み
	・画素値が0じゃない場所から光子を出す
		・ピクセル内でランダムの位置から
		・ランダムな角度で
	・検出（感度マップの完成）

	・

	・正規化(最大値で割る)
	・逆数をとってとりあえず画像表示
	・これを実際の投影画像に掛ける
	*/

	Condition cond;
	cond.detector_num = 1;
	cond.detector_size_w = 180;
	cond.detector_size_h = 180;
	cond.photon_num = 1000000;
	cond.surface_source_w = 2048;
	cond.surface_source_h = 1024;
	cond.distance_collimator_to_detector = 7.5;
	cond.pixel_size_detector = 0.2;
	cond.collimator_h = 1.;
	cond.collimator_w = 0.5;

	std::ostringstream ostr;

	float* surface_source = (float*)calloc(cond.surface_source_w * cond.surface_source_w, sizeof(float));
	float* efficiency_map = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));

	std::string surface_source_name;
	ostr << "./mksurface_source_img/surfaceSource_float_" << cond.surface_source_w <<  "-" << cond.surface_source_h << ".raw";
	surface_source_name = ostr.str();
	ostr.str("");
	readRawFile(surface_source_name, cond.surface_source_w * cond.surface_source_w, surface_source);

	mkEfficiencyMap(surface_source, efficiency_map, cond);


	std::string efficiency_map_name;
	ostr << "./result/efficiency_map_float_" << cond.detector_size_w <<  "-" << cond.detector_size_h << "-" << cond.detector_num << ".raw";
	efficiency_map_name = ostr.str();
	ostr.str("");
	writeRawFile(efficiency_map_name, cond.detector_num * cond.detector_size_w * cond.detector_size_h, efficiency_map);


	// float* sensitivity_map = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	// float* proj_img = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	//
	// readRawFile("primary_single_pinhole_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, sensitivity_map);
	// readRawFile("primary_single_pinhole_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, proj_img);
	// MakeSensitivityMap(sensitivity_map, cond);
	// SensitivityCorrection(sensitivity_map, proj_img, cond);
	//
	// writeRawFile("filter_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, proj_img);

}


void mkEfficiencyMap(float* surface_source, float* efficiency_map, struct Condition cond)
{
	init_genrand((unsigned)time(NULL));
	for(int theta_degree = 0; theta_degree < 360; theta_degree += 360 / cond.detector_num)
	{
		for(int i = 0; i < cond.surface_source_h; i++)
		{
			std::cout << "processing... " << i << " / " << cond.surface_source_h << std::endl;
			for(int j = 0; j < cond.surface_source_w; j++)
			{
				// 面線源画像の画素値が0の時は光子を飛ばさない
				if(surface_source[i * cond.surface_source_w + j] < 0.0001) { continue; }

				for(int p_num = 0; p_num < cond.photon_num; p_num++)
				{
					// 飛ぶ方向決定
					Photon p(i, j, cond);
					// if(abs(p.curr_(0)) > 0.4 && abs(p.curr_(2)) > 0.4) { continue; }
					// std::cout << "-- after setting photon -- \n";
					// std::cout << "x = " << p.curr_(0) << ", y = " << p.curr_(1) << ", z = " << p.curr_(2) << std::endl;

					Eigen::Vector3f on_detector;

					float abs_y = cond.distance_collimator_to_detector + cond.collimator_h / 2.;
					on_detector(0) = (p.curr_(0) + abs_y * tan(p.phi));
					on_detector(1) = - abs_y;
					on_detector(2) = (p.curr_(2) + abs_y * tan(M_PI / 2. - p.theta));

					// std::cout << "-- on detector -- \n";
					// std::cout << "x = " << on_detector(0) << ", y = " << on_detector(1) << ", z = " << on_detector(2) << std::endl;
					// std::cout << "\n";

					// コリメータの範囲内かどうかを検出
					int passible = canPassCollimator(p, on_detector, theta_degree, cond);
					if(!passible) { continue; }

					// 該当するコリメータで検出
					//検出器のピクセルサイズが0.5 cmであるため / 0.5  (×2)をしている
					// yの値が大きい→検出器の番号は小さい
					on_detector(0) /= cond.pixel_size_detector;
					on_detector(2) /= cond.pixel_size_detector;

					float j0 = (cond.detector_size_w - 1.) / 2. + on_detector(0);
					// float i0 = - cond.detector_size_h / 2. + on_detector(2);
					float i0 = (cond.detector_size_h - 1.) / 2. - on_detector(2);
					if(j0 < 0. || j0 > cond.detector_size_w) { continue; }
					if(i0 < 0. || i0 > cond.detector_size_h) { continue; }

					int j1 = (int)floor(j0);
					int i1 = (int)floor(i0);

					// positionとenergyを検出
					int index = (theta_degree / (360 / cond.detector_num)) * cond.detector_size_w * cond.detector_size_h + i1 * cond.detector_size_w + j1;
					// std::cout << index << std::endl;

					efficiency_map[index] += 1;
				}
			}
		}
	}
}


int canPassCollimator(class Photon p, Eigen::Vector3f on_detector, int theta_degree, struct Condition cond)
{
	Eigen::Vector3f on_colimator;

	// コリメータ手前
	on_colimator(1) = 0.;
	on_colimator(0) = p.curr_(0);
	on_colimator(2) = p.curr_(2);
	float colimator_radius = cond.collimator_w / 2. + tan(M_PI / 6);

	if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { return 0; }

	// コリメータ奥
	on_colimator(1) = - cond.collimator_h;
	float vec_scale = (on_colimator(1) - p.curr_(1)) / (on_detector(1) - p.curr_(1));
	on_colimator(0) = p.curr_(0) + vec_scale * (on_detector(0) - p.curr_(0));
	on_colimator(2) = p.curr_(2) + vec_scale * (on_detector(2) - p.curr_(2));
	colimator_radius = cond.collimator_w / 2. + tan(M_PI / 6);

	if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { return 0; }

	// コリメータ真ん中
	on_colimator(1) = - cond.collimator_h / 2.;
	vec_scale = (on_colimator(1) - p.curr_(1)) / (on_detector(1) - p.curr_(1));
	on_colimator(0) = p.curr_(0) + vec_scale * (on_detector(0) - p.curr_(0));
	on_colimator(2) = p.curr_(2) + vec_scale * (on_detector(2) - p.curr_(2));
	colimator_radius = cond.collimator_w / 2.;

	if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.)) { return 0; }
	// std::cout << "processing..." << std::endl;

	return 1;
}


void SensitivityCorrection(float* sensitivity_map, float* proj_img, struct Condition cond)
{
	// MakeSensitivityMap(sensitivity_map, cond);
	int pixel_num = cond.detector_num * cond.detector_size_w * cond.detector_size_h;

	for(int i = 0; i < pixel_num; i++) { proj_img[i] *= sensitivity_map[i]; }
}

void MakeSensitivityMap(float* sensitivity_map, struct Condition cond)
{
	int pixel_num = cond.detector_num * cond.detector_size_w * cond.detector_size_h;
	float max = 0.;
	for(int i = 0; i < pixel_num; i++) { if(max < sensitivity_map[i]) { max = sensitivity_map[i]; } }

	for(int i = 0; i < pixel_num; i++) { sensitivity_map[i] /= max; }

	for(int i = 0; i < pixel_num; i++) { if(sensitivity_map[i] > 0.01) { sensitivity_map[i] = 1. / sensitivity_map[i]; } }

	// 感度補正フィルタの確認
	std::ostringstream ostr;
	std::string efficiency_filter_name;
	ostr << "./result/sensitivity_filter_float_" << cond.detector_size_w <<  "-" << cond.detector_size_h << ".raw";
	efficiency_filter_name = ostr.str();
	ostr.str("");
	writeRawFile(efficiency_filter_name, pixel_num, sensitivity_map);

}


template <class T>
void readRawFile (std::string fname, const size_t num, T* image)
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
void writeRawFile(std::string fname, const size_t num, T* image)
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
