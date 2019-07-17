/*

floatにしてある
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)
*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "fileio.h"

typedef struct {
	int detector_size_w;
	int detector_size_h;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float d_width;
	float d_height;
	float time;
	float photon_num;
} Condition;

bool JudgeIsair(Eigen::Vector3f curr_, float* absorp_map, float img_pixel_size, int img_w, int img_h, int img_d);
float RandGenerator(curandState* state_gpu);
void readXcom(thrust::host_vector<float>& Coherent_ca, thrust::host_vector<float>& Compton_ca, thrust::host_vector<float>& Photoelectric_ca, thrust::host_vector<float>& mu_ca, thrust::host_vector<float>& Coherent_h2o, thrust::host_vector<float>& Compton_h2o, thrust::host_vector<float>& Photoelectric_h2o, thrust::host_vector<float>& mu_h2o);
void WriteImage(thrust::host_vector<float> &detector, Condition cond);
void WriteImage(thrust::host_vector<float> &detector, Condition cond, std::string directory);
void outputLog(Condition cond);
void outputLog(Condition cond, thrust::host_vector<float> &detector);
void make_directory(Condition cond, float w_collimator_mm);
std::string showCurrentTime();
std::string makeLogDirectory();

void launch_test_gradient_pinhole(std::vector<flaot> &detector, Condition cond);
void detect_photon(Photon p, std::vector<float> &detector, Condition cond);
bool passPinhole(Photon p, Condition cond);




class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float theta_;
	float phi_;
	Photon();
};

Photon::Photon()
{
	float rnd = 1. * 2 * (curand_uniform(&st) - 0.5f);
	theta_ = acos(rnd);
	phi_ = curand_uniform(&st) * 2 * M_PI;

	past_ << 0., 0., 0.;
	curr_ << 0., 0., 0.;
}

int main()
{
	Condition cond;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 13;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.2;
	cond.update_count = 10;
	cond.img_pixel_size = 0.2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.photon_num = 1000;

	/*---- 時間計測開始&条件表示 start ----*/
    cudaEvent_t start, end;

    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
    outputLog(cond);
	/*---- 時間計測開始&条件表示  end -----*/

	// 変数の定義
	std::vector<flaot> detector(cond.detector_size_w * cond.detector_size_h, 0.);

	// 関数呼び出し
	launch_test_gradient_pinhole(detector, cond);

	/*----- 時間計測処理 start-----*/
  cudaEventRecord(end, 0);
  cudaEventSynchronize(end);
  cudaEventElapsedTime(&cond.time, start, end);
  cond.time /= 1000.;
  outputLog(cond, detector);
  std::cout << "time = " << cond.time << "[s]\n"
					  <<	"     = " << cond.time/ 60 << "[min]\n"
					  << "     = " << cond.time / 3600 << "[h]" << std::endl;
	/*----- 時間計測処理 end-------*/
}


void launch_test_gradient_pinhole(std::vector<flaot> &detector, Condition cond)
{
	for(int i = 0; i < cond.photon_num; i++)
	{
		Photon p;

		bool pass_pinhole = passPinhole(p, cond);

		if(!pass_pinhole) break;

		detect_photon(p, detector, cond);
	}
}

bool passPinhole(Photon p, Condition cond)
{
	// l1の通過判定
	// l3の通過判定
	// l2の通過判定
}

void detect_photon(Photon p, std::vector<float> &detector, Condition cond)
{
	// 検出処理
}

void WriteImage(thrust::host_vector<float> &detector, Condition cond)
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

void WriteImage(thrust::host_vector<float> &detector, Condition cond, std::string directory)
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
	ostr << directory << "detector_single_pinhole_float_" << cond.detector_size_w << "-" << cond.detector_size_h << "-" << cond.detector_num << ".raw";
	write_proj_img_name = ostr.str();
	ostr.str("");
	writeRawFile(write_proj_img_name, proj_img_sum);
	/*-------------------------------------------*/

	/*----- primaryの検出結果書き込み -----*/
	std::string write_proj_img_name_primary;
	ostr << directory << "primary_single_pinhole_float_" << cond.detector_size_w << "-" << cond.detector_size_h << "-" << cond.detector_num << ".raw";
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

void outputLog(Condition cond, thrust::host_vector<float> &detector)
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

 		ostr << "\ntime = " << cond.time << "[s]\n"
			 	 << "     = " << cond.time / 60 << "[min]\n"
			   << "     = " << cond.time / 3600 << "[h]";

	str = ostr.str();


	std::string date_directory = makeLogDirectory();
	WriteImage(detector, cond, date_directory);
	std::string log_text = date_directory + "log.txt";
	std::ofstream outputfile(log_text.c_str());
	outputfile << str;
	outputfile.close();
}

void outputLog(Condition cond)
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

	str = ostr.str();

	std::cout << str << std::endl;
}
std::string showCurrentTime()
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y %B %d %A %H:%M:%S", date);
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
