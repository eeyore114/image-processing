
#include <stdio.h>
#include <unistd.h>
#include <iostream>
using namespace std;
#include <stdlib.h>

template <typename T> void readRawFile(std::string fname, std::vector<T>& image);
template <typename T> void writeRawFile(std::string fname, std::vector<T>& image);

// gnuplotを開いてグラフを描けるかテスト
void test_plot();
// 保存してあるファイルをグラフ化する場合
void test_plot_file_data();
// 配列の値をグラフ化する場合
void test_plot_array();

int main()
{
  test_plot();
  test_plot_file_data();
  test_plot_array();

  return 0;
}

void test_plot()
{
  FILE *gp;

  // -persist しないとcloseした時にグラフが消えちゃう
  gp = popen("gnuplot -persist", "w");
  fprintf(gp,"set terminal png\n");
  fprintf(gp,"set output \"test_plot.png\"\n");
  fprintf(gp, "set grid ;\n");
  fprintf(gp, "plot sin(x) w l\n");
  fflush(gp);
  pclose(gp);
}

void test_plot_file_data()
{
  FILE *gp, *data;
  std::string data_file;
  // データは必ずスペース区切りにする（カンマ区切りの場合は対応が必要）
  data_file="rec.csv";
  gp = popen("gnuplot -persist","w");
  fprintf(gp,"set terminal png\n");
  fprintf(gp,"set output \"test_plot_file_data.png\"\n");
  fprintf(gp, "set xrange [0:64]\n");
  fprintf(gp, "set yrange [0:70]\n");
  fprintf(gp, "plot \"%s\"  with lines\n",data_file.c_str());
  fflush(gp); // gpを吐き出す

  pclose(gp);
}

void test_plot_array()
{
  std::vector<float> img(64*64*64);
  readRawFile("ML-EM100_float_64-64-64.raw", img);
  std::vector<float> f(64, 0.0f);
  for(int i = 0; i < 64; i++) f[i] = img[31*64*64 + 31 * 64 + i];
  
  FILE *gp;
  gp=popen("gnuplot -persist","w");
  fprintf(gp,"set terminal png\n");
  fprintf(gp,"set output \"test_plot_array.png\"\n");
  fprintf(gp,"plot \"-\" w l\n");
  for(int i = 0; i < 64; i++) fprintf(gp,"%d %f\n", i, f[i]);
  fprintf(gp,"e\n");

  
  fflush(gp); // gpを吐き出す

  pclose(gp);
}

template <typename T> void writeRawFile(std::string fname, std::vector<T>& image)
{
    std::ofstream ofs(fname.c_str());
    if (!ofs) std::cerr << "Failed to open " << fname << std::endl;
    else { ofs.write(reinterpret_cast<char*>(&image[0]), sizeof(T) * image.size()); }
}


template <typename T> void readRawFile(std::string fname, std::vector<T>& image)
{
    std::ifstream ifs(fname.c_str());
    if (!ifs) std::cerr << "Failed to open " << fname << std::endl;
    else { ifs.read(reinterpret_cast<char*>(&image[0]), sizeof(T) * image.size()); }
}