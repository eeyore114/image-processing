#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "fileio.h"

void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h, std::vector<float> &original_img, struct Condition cond);
void make_cij(std::vector<float> &f, struct Condition cond);
void Projecion_SinglePinhole_3d(std::vector<float> &f, std::vector<float> &g, std::vector<float> &original_img, struct Condition cond, int is_inverse = 0);
void search_LessThan30_3d(std::vector<float> &f, std::vector<float> &g, struct Condition cond);
void outputTextFile(struct Condition cond, int is_pre = 0);
void saveTextFile(struct Condition cond, int is_before = 0);
std::string showCurrentTime();

template <class T>
void readRawFile (std::string fname, const size_t num, T* image);
template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);

struct Condition {
	int img_w;
	int img_h;
	int img_d;
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int update_count;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float img_pixel_size;
	float d_width;
	float d_height;
	float time;
	std::string text_name;
};

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
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.5;
	cond.update_count = 1;
	cond.img_pixel_size = 0.25;
	cond.d_width = 0.08;
	cond.d_height = 0.08;
	cond.text_name = "condition.txt";

	int is_inverse = 1;

	/*---- 時間計測開始&過去のログ表示 start ----*/
	clock_t start = clock();
	int is_pre = 1;
	outputTextFile(cond, is_pre);

	int is_before = 1;
	saveTextFile(cond, is_before);
	outputTextFile(cond);
	/*---- 時間計測開始&過去のログ表示  end -----*/


	// std::string readFileName1 = "backproj_float_128-128-128.raw";
	// std::string readFileName2 = "./read_img/efficiency_correction_float_512-256-180.raw";
	std::string readFileName2 = "./read_img/efficiency_correction_float_512-256-180.raw";
	std::string readFileName_original = "./read_img/mu-map_h2o_sphere12_float_128-128-128.raw";

	std::vector<float> f(cond.img_w * cond.img_h * cond.img_d, 0.);
	std::vector<float> g(cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> h(cond.img_w * cond.img_h * cond.img_d, 0.);
	std::vector<float> original_img(cond.img_w * cond.img_h * cond.img_d, 0.);
	readRawFile(readFileName_original, original_img);


	readRawFile(readFileName2, g);

	///////////////////////////////////////
	std::vector<float> test_f(f.size(), 0.);
	readRawFile(readFileName_original, test_f);
	std::vector<float> test_g(g.size(), 0.);
	std::vector<float> test_h(f.size(), 0.);
	Projecion_SinglePinhole_3d(test_f, test_g, original_img, cond);
	Projecion_SinglePinhole_3d(test_g, test_h, original_img, cond, is_inverse);
	writeRawFile("./result/test/non-absorp_float_512-256-180.raw", test_g);
	writeRawFile("./result/test/non-absorp_float_128-128-128.raw", test_h);


	///////////////////////////////////////

	/* 原画像だけの場合 */
	// readRawFile(readFileName_original, f);
	// Projecion_SinglePinhole_3d(f, g, original_img, cond);
	// writeRawFile("test_float_512-256-180.raw", g);
	// for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { f[i] = 0.; }
	/////////////////////
	// Projecion_SinglePinhole_3d(g, f, original_img, cond, is_inverse);
	// writeRawFile("./result/test/f_float_128-128-128.raw", f);

	std::ostringstream ostr;

	for(int i = 0; i < cond.update_count; i++)
	{
		std::cout << i + 1 << " times processing..." << std::endl;
		mlem(f, g, h, original_img, cond);
		/*---- 全部の回数保存する場合 ----*/
		std::string write_file_name;
		ostr << "result/ML-EM" <<  i + 1 << "_float_128-128-128.raw";
		write_file_name = ostr.str();
		ostr.str("");
		writeRawFile(write_file_name, h);
		/*-----------------------------*/
	}

	/*---- 最後だけ保存する場合 ----*/
	// std::string last;
	// ostr << "result/ML-EM" <<  update_count << "_float_128-128-128.raw";
	// last = ostr.str();
	// writeRawFile(last, h);
	/*---------------------------*/

	/*---- 時間計測終了&ログ書き込み start ----*/
	clock_t end = clock();

  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  saveTextFile(cond);
  std::cout << "time = " << cond.time << "[s]\n"
  					 <<	"     = " << cond.time / 60 << "[min]\n"
  					 << "     = " << cond.time / 3600 << "[h]" << std::endl;
	/*---- 時間計測終了&ログ書き込み  end  ----*/

}

void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h, std::vector<float> &original_img, struct Condition cond)
{
	std::vector<float> cij_proj(cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> cij(cond.img_w * cond.img_h * cond.img_d, 0.);
	std::vector<float> f_proj(cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> ratio_gf(cond.detector_num * cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> r_bproj(cond.img_w * cond.img_h * cond.img_d, 0.);
	int is_inverse = 1;

	make_cij(cij_proj, cond);
	Projecion_SinglePinhole_3d(cij_proj, cij, original_img, cond, is_inverse);
	Projecion_SinglePinhole_3d(f, f_proj, original_img, cond);
	writeRawFile("./result/test/f_proj_float_512-256-180.raw", f_proj);

	for(int i = 0; i < cond.detector_num * cond.detector_size_w * cond.detector_size_h; i++) { ratio_gf[i] =  g[i] / f_proj[i]; }

	Projecion_SinglePinhole_3d(ratio_gf, r_bproj, original_img, cond, is_inverse);

	for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { h[i] = (r_bproj[i] * f[i])  / cij[i]; }

	for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { f[i] = h[i]; }
}

void make_cij(std::vector<float> &f, struct Condition cond)
{
	for(int i = 0; i < cond.detector_num * cond.detector_size_w * cond.detector_size_h; i++) { f[i] = 1.; }
}

void Projecion_SinglePinhole_3d(std::vector<float> &f, std::vector<float> &g, std::vector<float> &original_img, struct Condition cond, int is_inverse)
{
	std::vector<float> theta_collimator_xy(cond.detector_size_w * cond.detector_size_h, 0.);
	std::vector<float> theta_collimator_zy(cond.detector_size_w * cond.detector_size_h, 0.);

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy, cond);
	int delta_detector = 360 / cond.detector_num;
	for(int theta_degree = 0; theta_degree < 360; theta_degree += delta_detector)
	{
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -1. * (cond.rotation_radius);

		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		float z = 0.;

		for(int m = 0; m < cond.detector_size_h; m++)
		{
			for(int n = 0; n < cond.detector_size_w; n++)
			{
				int count_h2o = 0;
				int count_ca = 0;

				//スケーリング処理
				float X = x / cond.img_pixel_size;
				float Y = y / cond.img_pixel_size;
				float Z = z / cond.img_pixel_size;

				if(theta_collimator_xy[m * cond.detector_size_w + n] != 0.0f || theta_collimator_zy[m * cond.detector_size_w + n] != 0.0f)
				{
					for (int sp = 0; sp < 300; sp++)
					{
						if(-(cond.img_w - 1.0f) / 2.0f < X && X < (cond.img_w - 1.0f) / 2.0f &&
						   -(cond.img_h - 1.0f) / 2.0f < Y && Y < (cond.img_h - 1.0f) / 2.0f &&
						   -(cond.img_d - 1.0f) / 2.0f < Z && Z < (cond.img_d - 1.0f) / 2.0f )
						{
							//双線形補完を行う左上の画素の座標(x0,y0)
							float x0 = floor(X - 0.5f) + 0.5f;
							float y0 = ceil(Y - 0.5f) + 0.5f;
							float z0 = ceil(Z - 0.5f) + 0.5f;

							//i,jに戻す
							float J = x0 + (cond.img_w - 1.0f) / 2.0f;
							float I = (cond.img_h - 1.0f) / 2.0f - y0;
							float D = (cond.img_d - 1.0f) / 2.0f - z0;

							int index1 = cond.img_w * cond.img_h * D + cond.img_w * I + J;
							int index2 = cond.img_w * cond.img_h * D + cond.img_w * I + (J + 1);
							int index3 = cond.img_w * cond.img_h * D + cond.img_w * (I + 1) + J;
							int index4 = cond.img_w * cond.img_h * D + cond.img_w * (I + 1) + (J + 1);
							int index5 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * I + J;
							int index6 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * I + (J + 1);
							int index7 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * (I + 1) + J;
							int index8 = cond.img_w * cond.img_h * (D + 1) + cond.img_w * (I + 1) + (J + 1);

							//(X , Y)から左上までの距離dx, dy, dz
							float dx = fabs(X - x0);
							float dy = fabs(Y - y0);
							float dz = fabs(Z - z0);

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
								g[index1] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V8;
								g[index2] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V7;
								g[index3] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V6;
								g[index4] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V5;
								g[index5] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V4;
								g[index6] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V3;
								g[index7] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V2;
								g[index8] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] * V1;
							}
							else
							{

								float val = f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 + f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;
								int j = round(X * cond.img_pixel_size + (cond.img_w - 1.0f) / 2.0f);
								int i = round((cond.img_h - 1.0f) / 2.0f - Y * cond.img_pixel_size);
								int d = round((cond.img_d - 1.0f) / 2.0f - Z * cond.img_pixel_size);
								if(original_img[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.26) { count_ca++; }
								if(original_img[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.14) { count_h2o++; }
								std::cout << "count_h2o = " << count_h2o << std::endl;

								float mu_h2o = 0.1538;
								float mu_ca = 0.2748;

								g[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] += val * expf(- mu_h2o * count_h2o * cond.img_pixel_size - mu_ca * count_ca * cond.img_pixel_size);
								// g[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] += val;
							}
						}
						X +=  cos(theta_collimator_zy[cond.detector_size_w * m + n]) * sin(theta - theta_collimator_xy[cond.detector_size_w * m + n]);
						Y +=  cos(theta_collimator_zy[cond.detector_size_w * m + n]) * cos(theta - theta_collimator_xy[cond.detector_size_w * m + n]);
						Z +=  -sin(theta_collimator_zy[cond.detector_size_w * m + n]);
					}
				}
			}
		}
	}
}

void search_LessThan30_3d(std::vector<float> &f, std::vector<float> &g, struct Condition cond)
{
	//確認用
	std::vector<float> h(cond.detector_size_w * cond.detector_size_h, 0.);

	for(int i = 0; i < cond.detector_size_h; i++)
	{
		for(int j = 0; j < cond.detector_size_w; j++)
		{
			float i0 = i * cond.d_width + cond.d_width / 2.;
			float j0 = j * cond.d_height + cond.d_height / 2.;
			Eigen::Vector3f on_detector;
			on_detector(0) = - cond.detector_size_w * cond.d_width / 2. + j0;
			on_detector(1) = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);
			on_detector(2) =   cond.detector_size_h * cond.d_height / 2. - i0;
			float theta_Ditector_xy = atanf(on_detector(0) / cond.distance_collimator_to_detector);
			float theta_Ditector_zy = atanf(on_detector(2) / cond.distance_collimator_to_detector);

			Eigen::Vector3f center_colimator;
			center_colimator << 0., -1. * cond.rotation_radius, 0.;

			// コリメータ奥
			Eigen::Vector3f on_colimator;
			on_colimator(1) = -1. * (cond.rotation_radius + cond.height_collimator / 2.);
			float vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
			on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
			on_colimator(2) = on_detector(2) + vec_scale * (center_colimator(2) - on_detector(2));
			float colimator_radius = cond.width_collimator / 2. + (cond.height_collimator / 2.) * tan(M_PI / 6);

			if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.))
			{
				f[i * cond.detector_size_w + j] = 0.0f;
				g[i * cond.detector_size_w + j] = 0.0f;
				h[i * cond.detector_size_w + j] = 0.;//確認用
				continue;
			}

			// コリメータ手前
			on_colimator(1) = -1. * (cond.rotation_radius - cond.height_collimator / 2.);
			vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
			on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
			on_colimator(2) = on_detector(2) + vec_scale * (center_colimator(2) - on_detector(2));
			colimator_radius = cond.width_collimator / 2. + (cond.height_collimator / 2.) * tan(M_PI / 6);

			if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.))
			{
				f[i * cond.detector_size_w + j] = 0.0f;
				g[i * cond.detector_size_w + j] = 0.0f;
				h[i * cond.detector_size_w + j] = 0.;//確認用
				continue;
			}

			f[i * cond.detector_size_w + j] = theta_Ditector_xy;
			g[i * cond.detector_size_w + j] = theta_Ditector_zy;
			h[i * cond.detector_size_w + j] = 100.0f;//確認用
		}
	}
	std::ostringstream ostr;
	std::string write_file_name;
	ostr << "test_within_30_float_" <<  cond.detector_size_w << "-" << cond.detector_size_h << ".raw";
	write_file_name = ostr.str();
	ostr.str("");
	writeRawFile(write_file_name, h);
}

void outputTextFile(struct Condition cond, int is_pre)
{
	//ファイルの読み込み
  std::ifstream fin( cond.text_name );

  if( !fin ) { std::cout << "failed to open condition file. \n"; }

  std::stringstream strstream;
  strstream << fin.rdbuf();
  fin.close();

  //ファイルの内容をstringに入れる
  std::string condition = strstream.str();

  //ファイルの内容を出力する
  if(is_pre) { std::cout << "----------- previous condition ----------\n"; }
  else 			 { std::cout << "--------------- condition ---------------\n"; }
  std::cout << condition << std::endl;
  std::cout << "-----------------------------------------\n\n";
}

void saveTextFile(struct Condition cond, int is_before)
{
	std::ofstream outputfile(cond.text_name);

  outputfile << "date : " << showCurrentTime() << "\n\n"
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
  					 << "img_pixel_size = " << cond.img_pixel_size;

  if(is_before == 0)
	{
		outputfile << "\ntime = " << cond.time << "[s]\n"
  					 	 << "     = " << cond.time / 60 << "[min]\n"
  					   << "     = " << cond.time / 3600 << "[h]";
	}
  outputfile.close();
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
