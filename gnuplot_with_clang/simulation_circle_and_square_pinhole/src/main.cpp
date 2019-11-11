/*

取得するデータ
・投影データ（円）
・再構成（円）
・投影データ（円）
・再構成（円）

投影データは円，矩形，原画像で重ねる
感度補正後と前両方
再構成も円，矩形，原画像で重ねる
*/

#include <stdio.h>
#include <unistd.h>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/Geometry"
#include "../include/fileio.h"
#include "../include/util.h"
#include "../include/util.cpp"

typedef struct {
	int img_w;
	int img_h;
	int img_d;
	int detector_w;
	int detector_h;
	int detector_num;
} Condition;

void plot_proj(Condition cond);
void plot_proj_after_efficiency_correction(Condition cond);
void plot_reconst_img(Condition cond);
 

int main()
{
	Condition cond;
	cond.img_w = cond.img_h = cond.img_d = 128;
	cond.detector_w = 512;
	cond.detector_h = 256;
	cond.detector_num = 180;
	plot_proj(cond);
	plot_proj_after_efficiency_correction(cond);
	plot_reconst_img(cond);
}

void plot_proj(Condition cond)
{
	// 拡張子を抜いた名前
	string write_profile_name = "projection_profile";
	vector<float> circle_proj(cond.detector_w * cond.detector_h * cond.detector_num);
	readRawFile("./img/primary_detector_float_512-256-180.raw", circle_proj);
	vector<float> square_proj(circle_proj.size());
	readRawFile("./img/primary_detector_rectangle_float_512-256-180.raw", square_proj);

	vector<float> circle_profile(cond.detector_w);
	vector<float> square_profile(circle_profile.size());
	rep(j, cond.detector_w)
	{
		int index = (cond.detector_num / 2 - 1) * cond.detector_w * cond.detector_h + (cond.detector_h / 2 - 1) * cond.detector_w + j;
		circle_profile[j] = circle_proj[index];
		square_profile[j] = square_proj[index];
	}


	// その配列をgnuplotで表示
	FILE *gp;
  gp=popen("gnuplot -persist","w");
	fprintf(gp,"set terminal eps\n");
	ostringstream ostr;
	ostr << "./result/" << write_profile_name << ".eps";
	std::string eps_file = ostr.str(); ostr.str("");
  fprintf(gp,"set output \"%s\"\n", eps_file.c_str());
	fprintf(gp,"set xrange [210:300]\n");
	fprintf(gp,"set yrange [0:3200]\n");
  fprintf(gp,"set xlabel \'width\'\n");
  fprintf(gp,"set ylabel \'count\'\n");
  fprintf(gp,"set xlabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set ylabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set tics font \"Arial, 18\"\n");
  fprintf(gp,"set xlabel offset 0,-1\n");
  fprintf(gp,"set ylabel offset -2,0\n");
  fprintf(gp,"set bmargin 6\n");
  fprintf(gp,"set lmargin 14\n");
  fprintf(gp,"set key horizontal\n");
  fprintf(gp,"set key font \"Arial, 18\"\n");
	
	fprintf(gp,"plot \"-\" w l title \"circle\", \"-\" w l dt (10,5) title \"square\"\n");
  rep(i, cond.detector_w) fprintf(gp,"%d %f\n", i, circle_profile[i]);
  fprintf(gp,"e\n");
  rep(i, cond.detector_w) fprintf(gp,"%d %f\n", i, square_profile[i]);
  fprintf(gp,"e\n");

  fflush(gp); // gpを吐き出す
  pclose(gp);

	// .png だとfontの設定ができないため，epsで保存し，pngで変換 
	ostr << "convert -density 100 " << eps_file << " ./result/" << write_profile_name << ".png";
	std::string convert_command = ostr.str(); ostr.str("");
	system(convert_command.c_str());
	ostr << "rm " << eps_file;
	std::string rm_command = ostr.str(); ostr.str("");
	system(rm_command.c_str());
}

void plot_proj_after_efficiency_correction(Condition cond)
{
	// 拡張子を抜いた名前
	string write_profile_name = "projection_profile_after_efficiency_correction";
	vector<float> circle_proj(cond.detector_w * cond.detector_h * cond.detector_num);
	readRawFile("./img/efficiency_correction_float_512-256-180.raw", circle_proj);
	vector<float> square_proj(circle_proj.size());
	readRawFile("./img/efficiency_correction_square_float_512-256-180.raw", square_proj);

	vector<float> circle_profile(cond.detector_w);
	vector<float> square_profile(circle_profile.size());
	rep(j, cond.detector_w)
	{
		int index = (cond.detector_num / 2 - 1) * cond.detector_w * cond.detector_h + (cond.detector_h / 2 - 1) * cond.detector_w + j;
		circle_profile[j] = circle_proj[index];
		square_profile[j] = square_proj[index];
	}

	// その配列をgnuplotで表示
	FILE *gp;
  gp=popen("gnuplot -persist","w");
	fprintf(gp,"set terminal eps\n");
	ostringstream ostr;
	ostr << "./result/" << write_profile_name << ".eps";
	std::string eps_file = ostr.str(); ostr.str("");
  fprintf(gp,"set output \"%s\"\n", eps_file.c_str());
	fprintf(gp,"set xrange [210:300]\n");
	fprintf(gp,"set yrange [0:3700]\n");
  fprintf(gp,"set xlabel \'width\'\n");
  fprintf(gp,"set ylabel \'count\'\n");
  fprintf(gp,"set xlabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set ylabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set tics font \"Arial, 18\"\n");
  fprintf(gp,"set xlabel offset 0,-1\n");
  fprintf(gp,"set ylabel offset -2,0\n");
  fprintf(gp,"set bmargin 6\n");
  fprintf(gp,"set lmargin 14\n");
  fprintf(gp,"set key horizontal\n");
  fprintf(gp,"set key font \"Arial, 18\"\n");
	
	fprintf(gp,"plot \"-\" w l title \"circle\", \"-\" w l dt (10,5) title \"square\"\n");
  rep(i, cond.detector_w) fprintf(gp,"%d %f\n", i, circle_profile[i]);
  fprintf(gp,"e\n");
  rep(i, cond.detector_w) fprintf(gp,"%d %f\n", i, square_profile[i]);
  fprintf(gp,"e\n");

  fflush(gp); // gpを吐き出す
  pclose(gp);

	// .png だとfontの設定ができないため，epsで保存し，pngで変換 
	ostr << "convert -density 100 " << eps_file << " ./result/" << write_profile_name << ".png";
	std::string convert_command = ostr.str(); ostr.str("");
	system(convert_command.c_str());
	ostr << "rm " << eps_file;
	std::string rm_command = ostr.str(); ostr.str("");
	system(rm_command.c_str());
}

void plot_reconst_img(Condition cond)
{
	// 拡張子を抜いた名前
	string write_profile_name = "ML-EM_100";
	vector<float> circle_reconst(cond.img_w * cond.img_h * cond.img_d);
	readRawFile("./img/ML-EM100_float_128-128-128.raw", circle_reconst);
	vector<float> square_reconst(circle_reconst.size());
	readRawFile("./img/ML-EM100_square_float_128-128-128.raw", square_reconst);

	vector<float> circle_profile(cond.img_w);
	vector<float> square_profile(circle_profile.size());
	rep(j, cond.img_w)
	{
		int index = (cond.img_d / 2 - 1) * cond.img_w * cond.img_h + (cond.img_h / 2 - 1) * cond.img_w + j;
		circle_profile[j] = circle_reconst[index];
		square_profile[j] = square_reconst[index];
	}

	// その配列をgnuplotで表示
	FILE *gp;
  gp=popen("gnuplot -persist","w");
	fprintf(gp,"set terminal eps\n");
	ostringstream ostr;
	ostr << "./result/" << write_profile_name << ".eps";
	std::string eps_file = ostr.str(); ostr.str("");
  fprintf(gp,"set output \"%s\"\n", eps_file.c_str());
	fprintf(gp,"set xrange [0:127]\n");
	fprintf(gp,"set yrange [0:100]\n");
  fprintf(gp,"set xlabel \'width\'\n");
  fprintf(gp,"set ylabel \'count\'\n");
  fprintf(gp,"set xlabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set ylabel font \"Times-Roman,30\"\n");
  fprintf(gp,"set tics font \"Arial, 18\"\n");
  fprintf(gp,"set xlabel offset 0,-1\n");
  fprintf(gp,"set ylabel offset -2,0\n");
  fprintf(gp,"set bmargin 6\n");
  fprintf(gp,"set lmargin 14\n");
  fprintf(gp,"set key horizontal\n");
  fprintf(gp,"set key font \"Arial, 18\"\n");
	
	fprintf(gp,"plot \"-\" w l title \"circle\", \"-\" w l dt (10,5) title \"square\"\n");
  rep(i, cond.img_w) fprintf(gp,"%d %f\n", i, circle_profile[i]);
  fprintf(gp,"e\n");
  rep(i, cond.img_w) fprintf(gp,"%d %f\n", i, square_profile[i]);
  fprintf(gp,"e\n");

  fflush(gp); // gpを吐き出す
  pclose(gp);

	// .png だとfontの設定ができないため，epsで保存し，pngで変換 
	ostr << "convert -density 100 " << eps_file << " ./result/" << write_profile_name << ".png";
	std::string convert_command = ostr.str(); ostr.str("");
	system(convert_command.c_str());
	ostr << "rm " << eps_file;
	std::string rm_command = ostr.str(); ostr.str("");
	system(rm_command.c_str());
}
