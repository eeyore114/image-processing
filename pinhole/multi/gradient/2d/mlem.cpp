/*

6/3 mlem関数の変数実装途中


変数名変更

writeImageは修正が必要

baseベクトルを変えれば良い

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
	int detector_num;
	int detector_size_w;
	int pinhole_count;
	int update_count;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float collimator_width;
	float img_pixel_size;
	float detector_pixel_size_w;
	float collimator_interval;
	float collimator_theta;
	float fov_theta;
	float time;
} Condition;

Eigen::Vector2f calculate_unit_vector(Eigen::Vector2f past, Eigen::Vector2f curr);
void create_fov(std::vector<int> &fov, Condition cond);
void launchProjection(std::vector<float> &init_img, std::vector<float> &detector, Condition cond);
void launchBackProjection(std::vector<float> &detector, std::vector<float> &reconstruct_img, Condition cond);
void launchMLEM(std::vector<float> &detector, std::vector<float> &reconstruct_img,  Condition cond);
void Projection(std::vector<float> &f, std::vector<float> &g, std::vector<int> &fov, Condition cond, bool is_inverse = false);
void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h,  Condition cond);
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
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.rotation_radius = 25;
	cond.pinhole_count = 3;
	cond.update_count = 100;
	cond.distance_collimator_to_detector = 7.5;
	cond.collimator_width = 0.2;
	cond.img_pixel_size = 0.2;
	cond.detector_pixel_size_w = 0.08;
	cond.collimator_interval = 9.;
	cond.collimator_theta = M_PI / 6.;
	cond.fov_theta = M_PI / 6.;

	clock_t start = clock();
	std::string date_directory;
	outputLogInit(cond, date_directory);

	std::vector<float> init_img(cond.img_w * cond.img_h, 0.);
	readRawFile("./read_img/Shepp_float_128-128.raw", init_img);
	std::vector<float> detector(cond.detector_size_w * cond.detector_num, 0.);
	launchProjection(init_img, detector, cond);
	std::vector<float> reconstruct_img(cond.img_w * cond.img_h, 0.);
	launchMLEM(detector, reconstruct_img, cond);


	clock_t end = clock();

  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
	outputLogLast(cond, date_directory, reconstruct_img);

  std::cout << "time = " << cond.time << "[s]\n"
  					 <<	"     = " << cond.time / 60 << "[min]\n"
  					 << "     = " << cond.time / 3600 << "[h]" << std::endl;
}

void launchMLEM(std::vector<float> &detector, std::vector<float> &reconstruct_img,  Condition cond)
{
	std::vector<float> init_img(cond.img_w * cond.img_h, 1.);

	std::ostringstream ostr;

	for(int i = 0; i < cond.update_count; i++)
	{
		std::cout << i + 1 << " times processing..." << std::endl;
		mlem(init_img, detector, reconstruct_img, cond);

		/*---- 全部の回数保存する場合 ----*/
		std::string write_file_name;
		ostr << "result/ML-EM" <<  i + 1 << "_float_" << cond.img_w <<  "-" << cond.img_h <<  ".raw";
		write_file_name = ostr.str();
		ostr.str("");
		if((i + 1) % 10 == 0) { writeRawFile(write_file_name, reconstruct_img); }
		/*-----------------------------*/
	}

	/*---- 最後だけ保存する場合 ----*/
	// std::string last;
	// ostr << "result/ML-EM" <<  cond.update_count << "_float_" << cond.img_w <<  "-" << cond.img_h <<  ".raw";
	// last = ostr.str();
	// writeRawFile(last, reconstruct_img);
	/*---------------------------*/
}

void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h,  Condition cond)
{
	std::vector<float> cij_proj(g.size(), 0.);
	std::vector<float> cij(f.size(), 0.);
	std::vector<float> f_proj(g.size(), 0.);
	std::vector<float> ratio_gf(g.size(), 0.);
	std::vector<float> ratio_gf_bproj(f.size(), 0.);

	for(int i = 0; i < cij_proj.size(); i++) { cij_proj[i] = 1.; }
	launchBackProjection(cij_proj, cij, cond);
	writeRawFile("./result/cij_float_128-128.raw", cij);

	launchProjection(f, f_proj, cond);
	writeRawFile("./result/f_proj_float_512-180.raw", f_proj);

	for(int i = 0; i < g.size(); i++) { ratio_gf[i] =  (g[i] < 0.0001 || f_proj[i] < 0.0001) ? 0 : g[i] / f_proj[i]; }

	launchBackProjection(ratio_gf, ratio_gf_bproj, cond);
	writeRawFile("./result/ratio_gf_float_512-180.raw", ratio_gf);

	// この下二つfor文使わないでかける？？
	for(int i = 0; i < h.size(); i++) { h[i] = (ratio_gf_bproj[i] * f[i])  / cij[i]; }

	for(int i = 0; i < h.size(); i++) { f[i] = h[i]; }
}


void launchProjection(std::vector<float> &init_img, std::vector<float> &detector, Condition cond)
{
	// fov求める
	std::vector<int> fov(cond.detector_size_w, -1);
	create_fov(fov, cond);
	std::vector<float> fov_float(fov.size(), -1.);
	for(int i = 0; i < fov.size(); i++) fov_float[i] = fov[i];

	writeRawFile("./result/fov_float_512-1.raw", fov_float);

	// 投影
	Projection(init_img, detector, fov, cond);

	// 画像書き込み
	writeRawFile("./result/detector_float_512-180.raw", detector);

}

void launchBackProjection(std::vector<float> &detector, std::vector<float> &reconstruct_img, Condition cond)
{
	std::vector<int> fov(cond.detector_size_w, -1);
	create_fov(fov, cond);
	bool is_inverse = true;
	Projection(detector, reconstruct_img, fov, cond, is_inverse);

	writeRawFile("./result/reconstruct_img_float_128-128.raw", reconstruct_img);
}

void Projection(std::vector<float> &f, std::vector<float> &g, std::vector<int> &fov, Condition cond, bool is_inverse)
{
	float pinhole_x[3] = { -9., 0., 9. };
	int delta_detector = 360 / cond.detector_num;

	for(int theta_degree = 0; theta_degree < 360; theta_degree += delta_detector)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for(int m = 0; m < cond.detector_size_w; m++)
		{
			if(fov[m] == -1) { continue; }
			for(int pinhole_num = 0; pinhole_num < cond.pinhole_count; pinhole_num++)
			{
				if(fov[m] != pinhole_num) { continue; }
				float collimator_x = pinhole_x[fov[m]];
				float collimator_y = - cond.rotation_radius;
				Eigen::Vector2f on_collimator;
				on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
				on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);

				float detector_x = (- (cond.detector_size_w - 1.) / 2. + m) * cond.detector_pixel_size_w;
				float detector_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);

				Eigen::Vector2f on_detector;
				on_detector(0) = detector_x * cosf(-theta) - detector_y * sinf(-theta);
				on_detector(1) = detector_x * sinf(-theta) + detector_y * cosf(-theta);

				Eigen::Vector2f d = calculate_unit_vector(on_detector, on_collimator);
				Eigen::Vector2f sp = on_detector / cond.img_pixel_size;

				for (int sample_point = 0; sample_point < 300; sample_point++)
				{
					if(-(cond.img_w - 1.0f) / 2.0f > sp(0) || sp(0) > (cond.img_w - 1.0f) / 2.0f || -(cond.img_h - 1.0f) / 2.0f > sp(1) || sp(1) > (cond.img_h - 1.0f) / 2.0f )
					{ sp += d; continue; }

					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(sp(0) - 0.5) + 0.5;
					float y0 = ceil(sp(1) - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 + (cond.img_w - 1.0)/2.0;
					float I = (cond.img_h - 1.0)/2.0 - y0;
					if(I + 1 == cond.img_h || J + 1 == cond.img_w) { continue; }
					//indexは配列の番号
					int index1 = cond.img_w * I + J;
					int index2 = cond.img_w * I + J + 1;
					int index3 = cond.img_w * (I + 1) + J;
					int index4 = cond.img_w * (I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(sp(0) - x0);
					float dy = fabs(sp(1) - y0);

					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下

					if(is_inverse)
					{
						if(isnan(g[index1]) ||  isnan(g[index2]) || isnan(g[index3]) || isnan(g[index4]))
						{

							// ここら辺でバグってる（IJをみてみる）
							std::cout << "index1 = " << index1 << std::endl;
							std::cout << "index2 = " << index2 << std::endl;
							std::cout << "index3 = " << index3 << std::endl;
							std::cout << "index4 = " << index4 << std::endl;
							std::cout << "I = " << I << std::endl;
							std::cout << "J = " << J << std::endl;
						}
						g[index1] += f[cond.detector_size_w * (theta_degree / delta_detector) + m] * S4;
						g[index2] += f[cond.detector_size_w * (theta_degree / delta_detector) + m] * S3;
						g[index3] += f[cond.detector_size_w * (theta_degree / delta_detector) + m] * S2;
						g[index4] += f[cond.detector_size_w * (theta_degree / delta_detector) + m] * S1;
					}
					else
					{
						float val = f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;
						g[cond.detector_size_w * (theta_degree / delta_detector) + m] += val;
					}
					sp += d;
				}
			}
		}
	}
}


void create_fov(std::vector<int> &fov, Condition cond)
{
	float pinhole_x[3] = { -9., 0., 9. };

	float fov_radius = abs(cond.distance_collimator_to_detector * tan(cond.collimator_theta));
	float y = - (cond.rotation_radius + cond.distance_collimator_to_detector);

	for(int pinhole = 0; pinhole < cond.pinhole_count; pinhole++)
	{
		Eigen::Vector2f on_collimator;
		on_collimator << pinhole_x[pinhole], - cond.rotation_radius;
		Eigen::Vector2f on_detector;
		on_detector << pinhole_x[pinhole], y;
		Eigen::Vector2f origin;
		origin << 0., 0.;
		Eigen::Vector2f base = calculate_unit_vector(origin, on_collimator);
		// これを直す
		Eigen::Vector2f max_vec;
		// 回転行列
		max_vec(0) = base(0) * cos(cond.fov_theta) - base(1) * sin(cond.fov_theta);
		max_vec(1) = base(0) * sin(cond.fov_theta) + base(1) * cos(cond.fov_theta);

		for(int i = 0; i < cond.detector_size_w; i++)
		{
			Eigen::Vector2f vec;
			float x = (- (cond.detector_size_w - 1.0) / 2.0 + i) * cond.detector_pixel_size_w;
			vec << x, y;
			vec = calculate_unit_vector(on_collimator, vec);
			if(max_vec.dot(base) > vec.dot(base)) { continue; }
			fov[i] = pinhole;
		}
	}
}

Eigen::Vector2f calculate_unit_vector(Eigen::Vector2f past, Eigen::Vector2f curr)
{
	Eigen::Vector2f d = curr - past;
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
			 << "detector_num = " << cond.detector_num << "\n"
			 << "detector_size_w = " << cond.detector_size_w << "\n"
			 << "rotation_radius = " << cond.rotation_radius << "\n"
			 << "distance_collimator_to_detector = " << cond.distance_collimator_to_detector << "\n"
			 << "collimator_width = " << cond.collimator_width << "\n"
			 << "detector_pixel_size_w = " << cond.detector_pixel_size_w << "\n"
			 << "img_pixel_size = " << cond.img_pixel_size << "\n"
			 << "pinhole_count = " << cond.pinhole_count << "\n"
			 << "collimator_interval = " << cond.collimator_interval << "\n"
			 << "collimator_theta = " << cond.collimator_theta << "\n"
			 << "update_count = " << cond.update_count << "\n"
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
	ostr << directory << "ML-EM" <<  cond.update_count << "_float_" << cond.img_w <<  "-" << cond.img_h << ".raw";
	result_name = ostr.str();
	ostr.str("");
	writeRawFile(result_name, result);
}
