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
void Projecion_SinglePinhole_3d(float* f,float* g, struct Condition cond);
void backProjection_Singlepinhole_3d(float* f, float* g, struct Condition cond);
void search_LessThan30_3d(float* f, float* g, struct Condition cond);

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
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
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

	int update_count = 10;

	// std::string readFileName1 = "backproj_float_128-128-128.raw";
	std::string readFileName2 = "./read_img/efficiency_correction_float_512-256-180.raw";
	// std::string readFileName2 = "./read_img/primary_float_180-180-180.raw";
	std::string readFileName_test = "./read_img/mu-map_h2o_sphere_float_128-128-128.raw";



	float* f = (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));
	float* g = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* h = (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));
	unsigned char writeFileName[50];

	readRawFile(readFileName2, cond.detector_num * cond.detector_size_w * cond.detector_size_h, g);
	/* 原画像だけの場合 */
	// readRawFile(readFileName_test, cond.img_w * cond.img_h * cond.img_d, f);
	// Projecion_SinglePinhole_3d(f, g, cond);
	// writeRawFile("test_float_512-256-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, g);
	// for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { f[i] = 0.; }
	/////////////////////

	backProjection_Singlepinhole_3d(g, f, cond);
	writeRawFile("test_float_128-128-128.raw", cond.img_w * cond.img_h * cond.img_d, f);

	std::string write_file_name;
	std::ostringstream ostr;

	for(int i = 0; i < update_count; i++)
	{
		std::cout << i + 1 << " times processing..." << std::endl;
		mlem(f, g, h, cond);
		// ostr << "result/ML-EM" <<  i + 1 << ".raw";
		// write_file_name = ostr.str();
		// ostr.str("");
		// writeRawFile(write_file_name, cond.img_w * cond.img_h * cond.img_d, h);
	}

	std::string last;
	ostr << "result/ML-EM" <<  update_count << "_float_128-128-128.raw";
	last = ostr.str();
	writeRawFile(last, cond.img_w * cond.img_h * cond.img_d, h);
}

void mlem(float* f, float* g, float* h, struct Condition cond)
{
	float* cij_proj = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* cij 		= (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));
	float* f_proj 	= (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* ratio_gf = (float*)calloc(cond.detector_num * cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* r_bproj 	= (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));

	make_cij(cij_proj, cond);
	backProjection_Singlepinhole_3d(cij_proj, cij, cond);
	writeRawFile("85_float_128-128-128.raw", cond.img_w * cond.img_h * cond.img_d, cij);
	Projecion_SinglePinhole_3d(f, f_proj, cond);
	writeRawFile("87_float_180-180-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, f_proj);


	for(int i = 0; i < cond.detector_num * cond.detector_size_w * cond.detector_size_h; i++) { ratio_gf[i] =  g[i] / f_proj[i]; }
	writeRawFile("91_float_180-180-180.raw", cond.detector_num * cond.detector_size_w * cond.detector_size_h, ratio_gf);

	backProjection_Singlepinhole_3d(ratio_gf, r_bproj, cond);
	writeRawFile("94_float_128-128-128.raw", cond.img_w * cond.img_h * cond.img_d, r_bproj);


	for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { h[i] = (r_bproj[i] * f[i])  / cij[i]; }
	// writeRawFile("98_float_128-128-128.raw", cond.img_w * cond.img_h * cond.img_d, h);

	for(int i = 0; i < cond.img_w * cond.img_h * cond.img_d; i++) { f[i] = h[i]; }
}

void make_cij(float* f, struct Condition cond)
{
	for(int i = 0; i < cond.detector_num * cond.detector_size_w * cond.detector_size_h; i++) { f[i] = 1.; }
}

void Projecion_SinglePinhole_3d(float* f,float* g, struct Condition cond)
{
	float pixel_size = 0.25;
	float* theta_collimator_xy = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* theta_collimator_zy = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy, cond);

	for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
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
				//スケーリング処理
				float X = x / pixel_size;
				float Y = y / pixel_size;
				float Z = z / pixel_size;

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

							g[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] += f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 +
							 									         f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;
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

void backProjection_Singlepinhole_3d(float* f, float* g, struct Condition cond)
{
	//ピクセルサイズは0.17 cm × 0.17 cm
	// float pixel_size = 0.17f;

	float pixel_size = 0.25;
	float* theta_collimator_xy = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* theta_collimator_zy = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy, cond);

	for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
	{
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -1. * (cond.rotation_radius);


		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		Eigen::Vector3f on_collimator;
		on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		on_collimator(2) = 0.;

		for(int m = 0; m < cond.detector_size_h; m++)
		{
			for(int n = 0; n < cond.detector_size_w; n++)
			{
				//スケーリング処理
				float X = on_collimator(0) / pixel_size;
				float Y = on_collimator(1) / pixel_size;
				float Z = on_collimator(2) / pixel_size;

				if(abs(theta_collimator_xy[m * cond.detector_size_w + n]) > 0.0001 || abs(theta_collimator_zy[m * cond.detector_size_w + n]) > 0.0001)
				{
					for (int sp = 0; sp < 800; sp++)
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

							g[index1] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V8;
							g[index2] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V7;
							g[index3] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V6;
							g[index4] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V5;
							g[index5] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V4;
							g[index6] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V3;
							g[index7] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V2;
							g[index8] += f[cond.detector_size_w * cond.detector_size_h * (theta_degree / 2) + cond.detector_size_w * m + n] * V1;
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


void search_LessThan30_3d(float* f, float* g, struct Condition cond)
{
	//確認用
	float* h = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));

	// 検出器の幅
	float d_width = 0.08;
	float d_height = 0.08;
	for(int i = 0; i < cond.detector_size_h; i++)
	{
		for(int j = 0; j < cond.detector_size_w; j++)
		{
			float i0 = i * d_width + d_width / 2.;
			float j0 = j * d_height + d_height / 2.;
			Eigen::Vector3f on_detector;
			on_detector(0) = - cond.detector_size_w * d_width / 2. + j0;
			on_detector(1) = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);
			on_detector(2) =   cond.detector_size_h * d_height / 2. - i0;
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
			float colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

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
			colimator_radius = cond.width_collimator / 2. + tan(M_PI / 6);

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
	std::string write_file_name = "test_within_30_float_180-180.raw";
	writeRawFile(write_file_name, cond.detector_size_w * cond.detector_size_h, h);
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
