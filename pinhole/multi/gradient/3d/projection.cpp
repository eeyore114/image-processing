/*

6/3 mlem関数の変数実装途中


変数名変更

writeImageは修正が必要

*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <vector>
#include "fileio.h"
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"

typedef struct {
	int img_w;
	int img_h;
	int img_d;
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int pinhole_count;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float collimator_width;
	float collimator_height;
	float img_pixel_size;
	float detector_pixel_size_w;
	float detector_pixel_size_h;
	float collimator_interval;
	float collimator_theta;
	float time;
} Condition;

Eigen::Vector3f calculate_unit_vector(Eigen::Vector3f past, Eigen::Vector3f curr);
void create_fov(std::vector<int> &fov, std::vector<float> &pinhole_x, std::vector<float> &pinhole_z, Condition cond);
void launchProjection(std::vector<float> &init_img, std::vector<float> &detector, Condition cond);
void launchBackProjection(std::vector<float> &detector, std::vector<float> &reconstruct_img, Condition cond);
void Projection(std::vector<float> &f, std::vector<float> &g, std::vector<int> &fov,  std::vector<float> &pinhole_x, std::vector<float> &pinhole_z, Condition cond, bool is_inverse = false);
void outputLogInit(Condition cond, std::string &date_directory);
void outputLogLast(Condition cond, std::string date_directory, std::vector<float> &result);
void WriteImage(std::vector<float> &detector, Condition cond, std::string directory);
std::string makeLogDirectory();
std::string showCurrentTime();


int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.img_d = 128;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 25;
	cond.pinhole_count = 8;
	cond.distance_collimator_to_detector = 7.5;
	cond.collimator_width = 0.2;
	cond.collimator_height = 0.2;
	cond.img_pixel_size = 0.2;
	cond.detector_pixel_size_w = 0.08;
	cond.detector_pixel_size_h = 0.08;
	cond.collimator_interval = 9.;
	cond.collimator_theta = M_PI / 6.;

	clock_t start = clock();
	std::string date_directory;
	outputLogInit(cond, date_directory);

	std::vector<float> init_img(cond.img_w * cond.img_h * cond.img_d, 0.);
	readRawFile("./read_img/Shepp_float_128-128-128.raw", init_img);
	std::vector<float> detector(cond.detector_size_w * cond.detector_size_h * cond.detector_num, 0.);
	launchProjection(init_img, detector, cond);

	std::vector<float> reconstruct_img(init_img.size(), 0.);

	launchBackProjection(detector, reconstruct_img, cond);


	clock_t end = clock();

  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
	outputLogLast(cond, date_directory, detector);

  std::cout << "time = " << cond.time << "[s]\n"
  					 <<	"     = " << cond.time / 60 << "[min]\n"
  					 << "     = " << cond.time / 3600 << "[h]" << std::endl;
}

void launchProjection(std::vector<float> &init_img, std::vector<float> &detector, Condition cond)
{
	std::vector<float> pinhole_x{ -15., 0., 15., -7.5, 7.5, -15.,  0., 15. };
	std::vector<float> pinhole_z{ 	 5., 5.,  5.,  0., 0.,  -5., -5., -5. };

	// fov求める
	std::vector<int> fov(cond.detector_size_w * cond.detector_size_h, -1);
	create_fov(fov, pinhole_x, pinhole_z, cond);
	std::vector<float> fov_float(cond.detector_size_w * cond.detector_size_h);
	for(int i = 0; i < fov.size(); i++) { fov_float[i] = fov[i];}

	writeRawFile("./result/fov_float_512-256.raw", fov_float);


	// 投影
	Projection(init_img, detector, fov, pinhole_x, pinhole_z, cond);

	// 画像書き込み
	writeRawFile("./result/detector_float_512-256-180.raw", detector);
}

void launchBackProjection(std::vector<float> &detector, std::vector<float> &reconstruct_img, Condition cond)
{
	std::vector<float> pinhole_x{ -15., 0., 15., -7.5, 7.5, -15.,  0., 15. };
	std::vector<float> pinhole_z{ 	 5., 5.,  5.,  0., 0.,  -5., -5., -5. };

	// fov求める
	std::vector<int> fov(cond.detector_size_w * cond.detector_size_h, -1);
	create_fov(fov, pinhole_x, pinhole_z, cond);

	// 投影
	bool is_inverse = true;
	Projection(detector, reconstruct_img, fov, pinhole_x, pinhole_z, cond, is_inverse);

	// 画像書き込み
	writeRawFile("./result/backproj_float_128-128-128.raw", reconstruct_img);
}


void Projection(std::vector<float> &f, std::vector<float> &g, std::vector<int> &fov,  std::vector<float> &pinhole_x, std::vector<float> &pinhole_z, Condition cond, bool is_inverse)
{
	int delta_detector = 360 / cond.detector_num;
	for(int theta_degree = 0; theta_degree < 360; theta_degree += delta_detector)
	{
		std::cout << theta_degree << " degree processing..." << std::endl;
		const float theta = theta_degree * M_PI / 180.0f;
		for(int m = 0; m < cond.detector_size_h; m++)
		{
			for(int n = 0; n < cond.detector_size_w; n++)
			{

				if(fov[m * cond.detector_size_w + n] == -1) { continue; }
				for(int pinhole_num = 0; pinhole_num < cond.pinhole_count; pinhole_num++)
				{
					if(fov[m * cond.detector_size_w + n] != pinhole_num) { continue; }
					float collimator_x = pinhole_x[fov[m * cond.detector_size_w + n]];
					float collimator_y = - cond.rotation_radius;
					float collimator_z = pinhole_z[pinhole_num];
					Eigen::Vector3f on_collimator;
					on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
					on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
					on_collimator(2) = collimator_z;

					float detector_x = (- (cond.detector_size_w - 1.) / 2. + n) * cond.detector_pixel_size_w;
					float detector_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);

					Eigen::Vector3f on_detector;
					on_detector(0) = detector_x * cosf(-theta) - detector_y * sinf(-theta);
					on_detector(1) = detector_x * sinf(-theta) + detector_y * cosf(-theta);
					on_detector(2) = ((cond.detector_size_h - 1.) / 2. - m) * cond.detector_pixel_size_h;

					Eigen::Vector3f d = calculate_unit_vector(on_detector, on_collimator);
					Eigen::Vector3f sp = on_detector / cond.img_pixel_size;

					for (int sample_point = 0; sample_point < 300; sample_point++)
					{
						if(-(cond.img_w - 1.0f) / 2.0f >sp(0) || sp(0) >(cond.img_w - 1.0f) / 2.0f ||
						   -(cond.img_h - 1.0f) / 2.0f >sp(1) || sp(1) >(cond.img_h - 1.0f) / 2.0f ||
						   -(cond.img_d - 1.0f) / 2.0f >sp(2) || sp(2) >(cond.img_d - 1.0f) / 2.0f )
						{ sp += d; continue; }

						//双線形補完を行う左上の画素の座標(x0,y0)
						float x0 = floor(sp(0) - 0.5) + 0.5;
						float y0 = ceil(sp(1) - 0.5) + 0.5;
						float z0 = ceil(sp(2) - 0.5) + 0.5;

						//i,jに戻す
						float J = x0 + (cond.img_w - 1.0)/2.0;
						float I = (cond.img_h - 1.0)/2.0 - y0;
						float D = (cond.img_d - 1.0f) / 2.0f - z0;
						if(I + 1 == cond.img_h || J + 1 == cond.img_w || D + 1 == cond.img_d) { continue; }
						//indexは配列の番号
						int index1 = cond.img_w * cond.img_h * D + cond.img_w * I + J;
						int index2 = cond.img_w * cond.img_h * D + cond.img_w * I + (J + 1);
						int index3 = cond.img_w * cond.img_h * D + cond.img_w * (I + 1) + J;
						int index4 = cond.img_w * cond.img_h * D + cond.img_w * (I + 1) + (J + 1);
						int index5 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * I + J;
						int index6 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * I + (J + 1);
						int index7 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * (I + 1) + J;
						int index8 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * (I + 1) + (J + 1);

						float dx = fabs(sp(0) - x0);
						float dy = fabs(sp(1) - y0);
						float dz = fabs(sp(2) - z0);

						float V1 = dx * dy * dz;
						float V2 = (1.0f - dx) * dy * dz;
						float V3 = dx * (1.0f - dy) * dz;
						float V4 = (1.0f - dx) * (1.0f - dy) * dz;
						float V5 = dx * dy * (1.0f - dz);
						float V6 = (1.0f - dx) * dy * (1.0f - dz);
						float V7 = dx * (1.0f - dy) * (1.0f - dz);
						float V8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);

						if(is_inverse)
						{
							int f_index = cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n;

							g[index1] += f[f_index] * V8;
							g[index2] += f[f_index] * V7;
							g[index3] += f[f_index] * V6;
							g[index4] += f[f_index] * V5;
							g[index5] += f[f_index] * V4;
							g[index6] += f[f_index] * V3;
							g[index7] += f[f_index] * V2;
							g[index8] += f[f_index] * V1;
						}
						else
						{
							float val = f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 + f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;

							g[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] += val;
						}
						sp += d;
					}
				}
			}
		}
	}
}

void create_fov(std::vector<int> &fov, std::vector<float> &pinhole_x, std::vector<float> &pinhole_z, Condition cond)
{

	float fov_radius = abs(cond.distance_collimator_to_detector * tan(cond.collimator_theta));
	float y = - (cond.rotation_radius + cond.distance_collimator_to_detector);

	for(int pinhole_num = 0; pinhole_num < cond.pinhole_count; pinhole_num++)
	{
		Eigen::Vector3f on_collimator;
		on_collimator << pinhole_x[pinhole_num], - cond.rotation_radius, pinhole_z[pinhole_num];
		Eigen::Vector3f on_detector;
		on_detector << pinhole_x[pinhole_num], y, pinhole_z[pinhole_num];
		Eigen::Vector3f base = calculate_unit_vector(on_collimator, on_detector);
		Eigen::Vector3f max_vec;
		// zコメントアウトしたら半径ちゃんとした大きさになった（ベクトルあんまりわかってない）
		max_vec << fov_radius + pinhole_x[pinhole_num], y, /*fov_radius + */pinhole_z[pinhole_num];
		Eigen::Vector3f max = calculate_unit_vector(on_collimator, max_vec);

		for(int i = 0; i < cond.detector_size_h; i++)
		{
			for(int j = 0; j < cond.detector_size_w; j++)
			{
				Eigen::Vector3f vec;
				float x = (- (cond.detector_size_w - 1.0) / 2.0 + j) * cond.detector_pixel_size_w;
				float z = ((cond.detector_size_h - 1.0) / 2.0 - i) * cond.detector_pixel_size_h;
				vec << x, y, z;
				vec = calculate_unit_vector(on_collimator, vec);
				if(max.dot(base) > vec.dot(base)) { continue; }
				fov[i * cond.detector_size_w + j] = pinhole_num;
			}
		}
	}
}

Eigen::Vector3f calculate_unit_vector(Eigen::Vector3f past, Eigen::Vector3f curr)
{
	Eigen::Vector3f d = curr - past;
	float d_norm = d.norm();
	d /= d_norm;
	return d;
}


void outputLogInit(Condition cond, std::string &date_directory)
{
	std::ostringstream ostr;
	ostr << "--------------- condition ---------------\n"
			 << "date : " << showCurrentTime() << "\n\n"
			 << "img_w = " << cond.img_w << "\n"
			 << "img_h = " << cond.img_h << "\n"
			 << "img_d = " << cond.img_d << "\n"
			 << "detector_num = " << cond.detector_num << "\n"
			 << "detector_size_w = " << cond.detector_size_w << "\n"
			 << "detector_size_h = " << cond.detector_size_h << "\n"
			 << "rotation_radius = " << cond.rotation_radius << "\n"
			 << "distance_collimator_to_detector = " << cond.distance_collimator_to_detector << "\n"
			 << "collimator_width = " << cond.collimator_width << "\n"
			 << "detector_pixel_size_w = " << cond.detector_pixel_size_w << "\n"
			 << "detector_pixel_size_h = " << cond.detector_pixel_size_h << "\n"
			 << "img_pixel_size = " << cond.img_pixel_size << "\n"
			 << "pinhole_count = " << cond.pinhole_count << "\n"
			 << "collimator_interval = " << cond.collimator_interval << "\n"
			 << "collimator_theta = " << cond.collimator_theta << "\n"
			 << "-----------------------------------------\n\n";

	std::string str = ostr.str();
	std::cout << str << std::endl;

 	date_directory = makeLogDirectory();
 	std::string log_text = date_directory + "log.txt";
 	std::ofstream outputfile(log_text.c_str());
 	outputfile << str;
 	outputfile.close();
}

void outputLogLast(Condition cond, std::string date_directory, std::vector<float> &result)
{
	std::ostringstream ostr;
	ostr << "\ntime = " << cond.time << "[s]\n"
			 << "     = " << cond.time / 60 << "[min]\n"
			 << "     = " << cond.time / 3600 << "[h]\n";

	std::string str = ostr.str();
	std::string log_text = date_directory + "log.txt";
	std::ofstream outputfile(log_text.c_str(), std::ios::app);
	outputfile << str;
	outputfile.close();

	WriteImage(result, cond, date_directory);
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


void WriteImage(std::vector<float> &result, Condition cond, std::string directory)
{
	std::ostringstream ostr;

	std::string result_name;
	ostr << directory << "projection_float_" << cond.detector_size_w << "-" << cond.detector_num << ".raw";
	result_name = ostr.str();
	ostr.str("");
	writeRawFile(result_name, result);
}
