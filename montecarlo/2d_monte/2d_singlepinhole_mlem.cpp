#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"

void mlem(float* f, float* g, float* h, struct Condition cond);
void make_cij(float* f, struct Condition cond);
void Projecion_SinglePinhole(float* f,float* g, struct Condition cond);
void backProjection_Singlepinhole(float* f, float* g, struct Condition cond);
void search_LessThan30(float* f, struct Condition cond);

template <class T>
void readRawFile (std::string fname, const size_t num, T* image);
template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);

struct Condition {
	int img_w;
	int img_h;
	int detector_num;
	int detector_size;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float pixel_size_img;
};

int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.detector_num = 180;
	cond.detector_size = 180;
	cond.rotation_radius = 15.;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.5;
	cond.pixel_size_img = 1;

	int update_count = 20;

	std::string readFileName2 = "cylinder_h2o_single_pinhole_float_180-180.raw";

	float* f = (float*)calloc(cond.img_w * cond.img_h, sizeof(float));
	float* g = (float*)calloc(cond.detector_num * cond.detector_size, sizeof(float));
	float* h = (float*)calloc(cond.img_w * cond.img_h, sizeof(float));
	unsigned char writeFileName[50];

	readRawFile(readFileName2, cond.detector_num * cond.detector_size, g);

	backProjection_Singlepinhole(g, f, cond);
	writeRawFile("after_b_proj_float_128-128.raw", cond.img_w * cond.img_h, f);

	std::string write_file_name;
	std::ostringstream ostr;

	for(int i = 0; i < update_count; i++)
	{
		std::cout << i + 1 << " times processing..." << std::endl;
		mlem(f, g, h, cond);
		ostr << "result/ML-EM" <<  i + 1 << "_float_" << cond.img_w << "-" << cond.img_h << ".raw";
		write_file_name = ostr.str();
		ostr.str("");
		// writeRawFile(write_file_name, cond.img_w * cond.img_h, h);
	}
	std::string last;
	ostr << "result/ML-EM" <<  update_count << "_float_" << cond.img_w << "-" << cond.img_h << ".raw";
	last = ostr.str();
	writeRawFile(last, cond.img_w * cond.img_h, h);
}

void mlem(float* f, float* g, float* h, struct Condition cond)
{
	float* cij_proj = (float*)calloc(cond.detector_num * cond.detector_size, sizeof(float));
	float* cij 		= (float*)calloc(cond.img_w * cond.img_h, sizeof(float));
	float* f_proj 	= (float*)calloc(cond.detector_num * cond.detector_size, sizeof(float));
	float* ratio_gf = (float*)calloc(cond.detector_num * cond.detector_size, sizeof(float));
	float* r_bproj 	= (float*)calloc(cond.img_w * cond.img_h, sizeof(float));

	make_cij(cij_proj, cond);
	backProjection_Singlepinhole(cij_proj, cij, cond);
	Projecion_SinglePinhole(f, f_proj, cond);

	for(int i = 0; i < cond.detector_num * cond.detector_size; i++) { ratio_gf[i] =  g[i] / f_proj[i]; }

	backProjection_Singlepinhole(ratio_gf, r_bproj, cond);


	for(int i = 0; i < cond.img_w * cond.img_h; i++) { h[i] = (r_bproj[i] * f[i])  / cij[i]; }

	for(int i = 0; i < cond.img_w * cond.img_h; i++) { f[i] = h[i]; }
}

void make_cij(float* f, struct Condition cond)
{
	for(int i = 0; i < cond.detector_num * cond.detector_size; i++) { f[i] = 1.; }
}

void Projecion_SinglePinhole(float* f,float* g, struct Condition cond)
{
	float* theta_collimator = (float*)calloc(cond.detector_size, sizeof(float));

	search_LessThan30(theta_collimator, cond);

	for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
	{
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);


		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		Eigen::Vector2f on_collimator;
		on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);

		for(int m = 0; m < cond.detector_size; m++)
		{
			//スケーリング処理
			float X = on_collimator(0) / cond.pixel_size_img;
			float Y = on_collimator(1) / cond.pixel_size_img;

			if(theta_collimator[m] != 0.0f)
			{
				for (int sp = 0; sp < 300; sp++)
				{
					if(-(cond.img_w - 1.0f) / 2.0f < X && X < (cond.img_w - 1.0f) / 2.0f &&
					   -(cond.img_h - 1.0f) / 2.0f < Y && Y < (cond.img_h - 1.0f) / 2.0f  )
					{
						//双線形補完を行う左上の画素の座標(x0,y0)
						float x0 = floor(X - 0.5f) + 0.5f;
						float y0 = ceil(Y - 0.5f) + 0.5f;

						//i,jに戻す
						float J = x0 + (cond.img_w - 1.0f) / 2.0f;
						float I = (cond.img_h - 1.0f) / 2.0f - y0;

						//indexは配列の番号
						int index1 = cond.img_w * I + J;
						int index2 = cond.img_w * I + J + 1;
						int index3 = cond.img_w * (I + 1) + J;
						int index4 = cond.img_w * (I + 1) + J + 1;

						//(X , Y)から左上までの距離dx, dy, dz
						float dx = fabs(X - x0);
						float dy = fabs(Y - y0);

						//双線形補完で使う面積S1,S2,S3,S4
						float S1 = dx * dy;//左上
						float S2 = (1.0f - dx) * dy;//右上
						float S3 = dx * (1.0f - dy);//左下
						float S4 = (1.0f - dx) * (1.0f - dy);//右下

						
						g[cond.detector_size * theta_degree / (360 / cond.detector_num) + m] += f[index1] * S4 + f[index2] * S3 + f[index3] * S2 + f[index4] * S1;	
					}
					X +=  sin(theta - theta_collimator[m]);
					Y +=  cos(theta - theta_collimator[m]); 
				}
			}
		}
	}
}

void backProjection_Singlepinhole(float* f, float* g, struct Condition cond)
{

	float* theta_collimator = (float*)calloc(cond.detector_size, sizeof(float));

	search_LessThan30(theta_collimator, cond);

	for(int theta_degree = 0; theta_degree < 360; theta_degree += 360 / cond.detector_num)
	{
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);


		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		Eigen::Vector3f on_collimator;
		on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		on_collimator(2) = 0.;

		for(int m = 0; m < cond.detector_size; m++)
		{
			//スケーリング処理
			float X = on_collimator(0) / cond.pixel_size_img;
			float Y = on_collimator(1) / cond.pixel_size_img;

			if(theta_collimator[m] != 0.0f)
			{
				for (int sp = 0; sp < 300; sp++)
				{
					if(-(cond.img_w - 1.0f) / 2.0f < X && X < (cond.img_w - 1.0f) / 2.0f &&
					   -(cond.img_h - 1.0f) / 2.0f < Y && Y < (cond.img_h - 1.0f) / 2.0f  )
					{
						//双線形補完を行う左上の画素の座標(x0,y0)
						float x0 = floor(X - 0.5f) + 0.5f;
						float y0 = ceil(Y - 0.5f) + 0.5f;

						//i,jに戻す
						float J = x0 + (cond.img_w - 1.0f) / 2.0f;
						float I = (cond.img_h - 1.0f) / 2.0f - y0;

						//indexは配列の番号
						int index1 = cond.img_w * I + J;
						int index2 = cond.img_w * I + J + 1;
						int index3 = cond.img_w * (I + 1) + J;
						int index4 = cond.img_w * (I + 1) + J + 1;

						//(X , Y)から左上までの距離dx, dy, dz
						float dx = fabs(X - x0);
						float dy = fabs(Y - y0);

						//双線形補完で使う面積S1,S2,S3,S4
						float S1 = dx * dy;//左上
						float S2 = (1.0f - dx) * dy;//右上
						float S3 = dx * (1.0f - dy);//左下
						float S4 = (1.0f - dx) * (1.0f - dy);//右下

						g[index1] += f[cond.detector_size * theta_degree / (360 / cond.detector_num) + m] * S4;
						g[index2] += f[cond.detector_size * theta_degree / (360 / cond.detector_num) + m] * S3;
						g[index3] += f[cond.detector_size * theta_degree / (360 / cond.detector_num) + m] * S2;
						g[index4] += f[cond.detector_size * theta_degree / (360 / cond.detector_num) + m] * S1;
					}

					X +=  sin(theta - theta_collimator[m]);
					Y +=  cos(theta - theta_collimator[m]); 
				}
			}
		}
	}
}


void search_LessThan30(float* f, struct Condition cond)
{
	//確認用
	float* h = (float*)calloc(cond.detector_size, sizeof(float));

	// 検出器の幅
	float d_width = 0.2;
	for(int j = 0; j < cond.detector_size; j++)
	{
		
		float j0 = j * d_width + d_width / 2.;
		Eigen::Vector2f on_detector;
		on_detector(0) = - cond.detector_size * d_width / 2. + j0;
		on_detector(1) = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);
		float theta_Ditector_xy = atanf(on_detector(0) / cond.distance_collimator_to_detector);

		Eigen::Vector3f center_colimator;
		center_colimator << 0., -1. * cond.rotation_radius, 0.;

		// コリメータ奥
		Eigen::Vector2f on_colimator;
		on_colimator(1) = -1. * (cond.rotation_radius + cond.height_collimator / 2.);
		float vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
		on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
		float colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

		if(abs(on_colimator(0)) > colimator_radius)
		{
			f[j] = 0.0f;
			h[j] = 0.;//確認用
			continue;
		}

		// コリメータ手前
		on_colimator(1) = -1. * (cond.rotation_radius - cond.height_collimator / 2.);
		vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
		on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
		colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

		if(abs(on_colimator(0)) > colimator_radius)
		{
			f[j] = 0.0f;
			h[j] = 0.;//確認用
			continue;
		}

		f[j] = theta_Ditector_xy;
		h[j] = 100.0f;//確認用
	}
	std::string write_file_name = "test_within_30_float_65-1.raw";
	writeRawFile(write_file_name, cond.detector_size, h);
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
