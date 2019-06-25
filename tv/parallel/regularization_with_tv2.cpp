#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <vector>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "fileio.h"

void ParallelProjecion(std::vector<float> &f, std::vector<float> &g, struct Condition cond, int is_inverse = 0);
void RegurarizationWithTV(std::vector<float> &proj_img, std::vector<float> &bproj_img, struct Condition cond);
void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h, struct Condition cond);
void computeTVgradient(std::vector<float> &f, std::vector<float> &v, struct Condition cond);
void outputTextFile(struct Condition cond, int is_pre = 0);
void saveTextFile(struct Condition cond, int is_before = 0);

struct Condition {
	int img_w;
	int img_h;
	int detector_num;
	int detector_size;
	int update_count;
	float height_collimator;
	float width_collimator;
	float img_pixel_size;
	float d_width;
	float d_height;
	float distance_img_bottom_to_detector;
	float time;
	std::string text_name;
};

int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.detector_num = 20;
	cond.detector_size = 128;
	cond.update_count = 400;
	cond.distance_img_bottom_to_detector = 10;
	cond.text_name = "condition.txt";
	/*---- 使ってない ----*/
	cond.height_collimator = 1.;
	cond.width_collimator = 1.;
	cond.img_pixel_size = 1.;
	cond.d_width = 1.;
	cond.d_height = 1.;
	/*------------------*/

	clock_t start = clock();
	int is_pre = 1;
	outputTextFile(cond, is_pre);

	int is_before = 1;
	saveTextFile(cond, is_before);
	outputTextFile(cond);

	////////////////////////////////////////////////

	std::string read_original_img_name = "./read_img/shepp_float_128-128.raw";

	std::vector<float> original_img(cond.img_w * cond.img_h, 0.);
	std::vector<float> proj_img(cond.detector_num * cond.detector_size, 0.);
	std::vector<float> bproj_img(original_img.size(), 0.);
	readRawFile(read_original_img_name, original_img);

	ParallelProjecion(original_img, proj_img, cond);

	int is_inverse = 1;
	ParallelProjecion(proj_img, bproj_img, cond, is_inverse);

	RegurarizationWithTV(proj_img, bproj_img, cond);

	/*---- MLEMだけの場合 ----*/
	// std::vector<float> proj_img_only_mlem(cond.detector_num * cond.detector_size, 0.);
	// std::vector<float> bproj_img_only_mlem(original_img.size(), 0.);
	// std::vector<float> only_mlem_img(original_img.size(), 0.);
	// ParallelProjecion(original_img, proj_img_only_mlem, cond);
	// ParallelProjecion(proj_img_only_mlem, bproj_img_only_mlem, cond, is_inverse);
	// for(int i = 0; i < cond.update_count; i++) { mlem(bproj_img_only_mlem, proj_img_only_mlem, only_mlem_img, cond); }
	// std::ostringstream ostr;
	// std::string write_only_mlem_file_name;
	// ostr << "result/only_mlem" <<  cond.update_count << "_float_" << cond.img_w << "-" << cond.img_h << ".raw";
	// write_only_mlem_file_name = ostr.str();
	// ostr.str("");
	// writeRawFile(write_only_mlem_file_name, only_mlem_img);
	/*----------------------*/

	clock_t end = clock();

  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  saveTextFile(cond);
  std::cout << "time = " << cond.time << "[s]\n"
  					 <<	"     = " << cond.time / 60 << "[min]\n"
  					 << "     = " << cond.time / 3600 << "[h]" << std::endl;
}

void RegurarizationWithTV(std::vector<float> &proj_img, std::vector<float> &bproj_img, struct Condition cond)
{
	std::ostringstream ostr;
	int mlem_count = 5;
	int tv_count = 2;

	for(int update_count = 0; update_count < cond.update_count; update_count++)
	{
		std::cout << update_count + 1 << " times processing..." << std::endl;
		std::vector<float> bproj_img_with_tv(bproj_img.size(), 0.);
		for(int i = 0; i < mlem_count; i++) {	mlem(bproj_img, proj_img, bproj_img_with_tv, cond); }

		for(int i = 0; i < bproj_img_with_tv.size(); i++) { if(bproj_img_with_tv[i] < 0.) { bproj_img_with_tv[i] = 0.; }}

		std::vector<float> v(bproj_img.size(), 0.);
		for(int i = 0; i < tv_count; i++)
		{
			computeTVgradient(bproj_img_with_tv, v, cond);
			float a = 0.2;
			for(int i = 0; i < bproj_img.size(); i++) { bproj_img_with_tv[i] -= a * v[i]; }
		}
		for(int i = 0; i < bproj_img_with_tv.size(); i++) { if(bproj_img_with_tv[i] < 0.) { bproj_img_with_tv[i] = 0.; }}

		bproj_img = bproj_img_with_tv;


		/*---- 全部の回数保存する場合 ----*/
		// std::string write_file_name;
		// ostr << "result/with_tv" <<  update_count + 1 << "_float_" << cond.img_w << "-" << cond.img_h << ".raw";
		// write_file_name = ostr.str();
		// ostr.str("");
		// writeRawFile(write_file_name, bproj_img);
		/*-----------------------------*/
	}
	/*---- 全部の回数保存する場合 ----*/
	std::string write_file_name;
	ostr << "result/with_tv" <<  cond.update_count << "_float_" << cond.img_w << "-" << cond.img_h << ".raw";
	write_file_name = ostr.str();
	ostr.str("");
	writeRawFile(write_file_name, bproj_img);
	/*-----------------------------*/
}

void ParallelProjecion(std::vector<float> &f, std::vector<float> &g, struct Condition cond, int is_inverse)
{
	int detect_sub = 360 / cond.detector_num;
	for(int theta_degree = 0; theta_degree < 360; theta_degree += detect_sub)
	{
		const float theta = theta_degree * M_PI / 180.0f;
		for (int i = 0; i < cond.img_h + cond.distance_img_bottom_to_detector + 10; i++)
		{
			for (int j = 0; j < cond.img_w; j++)
			{
				//検出器の座標を定義
				float x = - cond.img_w / 2 + j;
				float y = - cond.img_h / 2 - cond.distance_img_bottom_to_detector + i;

				//座標(x, y)を原点を中心に-theta回転させた座標(X, Y)
				float X = x * cosf(-theta) - y * sinf(-theta);
				float Y = x * sinf(-theta) + y * cosf(-theta);

				if(-(cond.img_w - 1.0f) / 2.0f < X && X < (cond.img_w - 1.0f) / 2.0f && -(cond.img_h - 1.0f) / 2.0f < Y && Y < (cond.img_h - 1.0f) / 2.0f )
				{
					//双線形補完を行う左上の画素の座標(x0,y0)
					float x0 = floor(X - 0.5) + 0.5;
					float y0 = ceil(Y - 0.5) + 0.5;

					//i,jに戻す
					float J = x0 +(cond.img_w - 1.0)/2.0 ;
					float I = (cond.img_h - 1.0)/2.0 - y0;

					//indexは配列の番号
					int index1 = cond.img_w * I + J;
					int index2 = cond.img_w * I + J + 1;
					int index3 = cond.img_w * ( I + 1) + J;
					int index4 = cond.img_w * ( I + 1) + J + 1;

					//(X , Y)から左上までの距離dx,dy
					float dx = fabs(X - x0);
					float dy = fabs(Y - y0);

					//双線形補完で使う面積S1,S2,S3,S4
					float S1 = dx * dy;//左上
					float S2 = (1.0f - dx) * dy;//右上
					float S3 = dx * (1.0f - dy);//左下
					float S4 = (1.0f - dx) * (1.0f - dy);//右下
					if(is_inverse == 0)
					{
						g[cond.img_w * (theta_degree / detect_sub) + j] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;
					}
					else
					{
						g[index1] += S4 * f[cond.img_w * (theta_degree / detect_sub) + j];
						g[index2] += S3 * f[cond.img_w * (theta_degree / detect_sub) + j];
						g[index3] += S2 * f[cond.img_w * (theta_degree / detect_sub) + j];
						g[index4] += S1 * f[cond.img_w * (theta_degree / detect_sub) + j];
					}
				}
			}
		}
	}
}

void mlem(std::vector<float> &f, std::vector<float> &g, std::vector<float> &h, struct Condition cond)
{
	std::vector<float> cij_proj(g.size(), 0.);
	std::vector<float> cij(f.size(), 0.);
	std::vector<float> f_proj(g.size(), 0.);
	std::vector<float> ratio_gf(g.size(), 0.);
	std::vector<float> r_bproj(f.size(), 0.);

	int is_inverse = 1;
	for(int i = 0; i < cij_proj.size(); i++) { cij_proj[i] = 1.; }
	ParallelProjecion(cij_proj, cij, cond, is_inverse);

	ParallelProjecion(f, f_proj, cond);

	for(int i = 0; i < ratio_gf.size(); i++) { ratio_gf[i] =  (g[i] < 0.0001 || f_proj[i] < 0.0001) ? 0 : g[i] / f_proj[i]; }

	ParallelProjecion(ratio_gf, r_bproj, cond, is_inverse);

	for(int i = 0; i < h.size(); i++) { h[i] = (r_bproj[i] * f[i])  / cij[i]; }

	for(int i = 0; i < f.size(); i++) {	f[i] = h[i]; }
}


void computeTVgradient(std::vector<float> &f, std::vector<float> &v, struct Condition cond)
{
	float epc = 0.0001;
	const float W = cond.img_w;
	for (int i = 0; i < cond.img_h; i++)
	{
		for (int j = 0; j < cond.img_w; j++)
		{
			if(i == 0 || i == cond.img_h - 1 || j == 0 || j == cond.img_w) { continue; }

			float term1 = ( 2. * (f[i * W + j] - f[(i - 1) * W + j]) + 2. * (f[i * W + j] - f[i * W + (j - 1)]) ) / sqrt(epc + pow(f[i * W + j] - f[(i - 1) * W + j], 2.) + pow(f[i * W + j] - f[i * W + (j - 1)], 2.));
			float term2 = ( 2. * (f[(i + 1) * W + j] - f[i * W + j]) ) / sqrt(epc + pow(f[(i + 1) * W + j] - f[i * W + j], 2.) + pow(f[(i + 1) * W + j] - f[(i + 1) * W + (j - 1)], 2.));
			float term3 = ( 2. * (f[i * W + (j + 1)] - f[i * W + j]) ) / sqrt(epc + pow(f[i * W + (j + 1)] - f[i * W + j], 2.) + pow(f[i * W + (j + 1)] - f[(i - 1) * W + (j + 1)], 2.));
			v[i * cond.img_w + j] = term1 - term2 - term3;
			// std::cout << "term1" << term1 << std::endl;
			// std::cout << "term2" << term2 << std::endl;
			// std::cout << "term3" << term3 << std::endl;
		}
	}
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
  if(is_pre) { std::cout << "--- previous condition ---\n"; }
  else 			 { std::cout << "------- condition -------\n"; }
  std::cout << condition << std::endl;
  std::cout << "-------------------------\n\n";
}

void saveTextFile(struct Condition cond, int is_before)
{
	std::ofstream outputfile(cond.text_name);

  outputfile << "img_w = " << cond.img_w << "\n"
  					 << "img_h = " << cond.img_h << "\n"
  					 << "detector_num = " << cond.detector_num << "\n"
  					 << "detector_size = " << cond.detector_size << "\n"
  					 << "update_count = " << cond.update_count << "\n"
  					 << "distance_img_bottom_to_detector = " << cond.distance_img_bottom_to_detector;

  if(is_before == 0)
	{
		outputfile << "\ntime = " << cond.time << "[s]\n"
  					 	 << "     = " << cond.time / 60 << "[min]\n"
  					   << "     = " << cond.time / 3600 << "[h]";
	}
  outputfile.close();
}
