#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>


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
	cond.update_count = 10;
	cond.text_name = "condition.txt";


	clock_t start = clock();


  //ファイルの読み込み
  std::ifstream fin( cond.text_name );

  if( !fin ) { std::cout << "failed to open condition file. \n"; }

  std::stringstream strstream;
  strstream << fin.rdbuf();
  fin.close();

  //ファイルの内容をstringに入れる
  std::string data = strstream.str();

  //ファイルの内容を出力する
  std::cout << data << std::endl;

  std::ofstream outputfile(cond.text_name);
  outputfile << "img_w = " << cond.img_w << "\n"
  					 << "img_h = " << cond.img_h << "\n"
  					 << "img_d = " << cond.img_d << "\n"
  					 << "detector_num = " << cond.detector_num << "\n"
  					 << "detector_size_w = " << cond.detector_size_w << "\n"
  					 << "rotation_radius = " << cond.rotation_radius << "\n";
  outputfile.close();

  clock_t end = clock();

  cond.time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  printf("time %lf[s]\n", cond.time);

  return 0;
}
