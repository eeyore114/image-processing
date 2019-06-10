/*

floatにしてある

genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

円柱ファントムの拡大比を
float ratio
で設定


*/

#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>

void SensitivityCorrection(float* f, float* g, struct Condition cond);
void MakeSensitivityMap(float* f, struct Condition cond);

struct Condition {
	int detector_num;
	int detector_size_w;
	int detector_size_h;
};


template <class T>
void readRawFile (string fname, const size_t num, T* image);

template <class T>
void writeRawFile (string fname, const size_t num, T* image);

int main()
{
	/*
	・点線源から出した投影画像を読み込み
	・正規化(最大値で割る)
	・逆数をとってとりあえず画像表示
	・これを実際の投影画像に掛ける
	*/

	Condition cond;
	cond.detector_num = 180;
	cond.detector_size_w = 50;
	cond.detector_size_h = 50;


	float* sensitivity_map = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* proj_img = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));

	readRawFile("primary_single_pinhole_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, sensitivity_map);
	readRawFile("primary_single_pinhole_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, proj_img);
	MakeSensitivityMap(sensitivity_map, cond);
	SensitivityCorrection(sensitivity_map, proj_img, cond);

	writeRawFile("filter_float_50-50-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, proj_img);

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

	// 感度補正フィルタの確認
	writeRawFile("sensitivity_filter_float_50-50-180.raw", pixel_num, sensitivity_map);

	for(int i = 0; i < pixel_num; i++) { if(sensitivity_map[i] > 0.01) { sensitivity_map[i] = 1. / sensitivity_map[i]; } }

	writeRawFile("inverse_sensitivity_filter_float_50-50-180.raw", pixel_num, sensitivity_map);
	
}


template <class T>
void readRawFile (string fname, const size_t num, T* image)
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
void writeRawFile(string fname, const size_t num, T* image)
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
