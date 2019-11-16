#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
	int update_count;
	int thread_num;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float img_pixel_size;
	float d_width;
	float d_height;
	float time;
} Condition;

const char readFileName[] = "Shepp_float_64-64-64.raw";
const char writeFileName1[] = "sphere_singlepinhole_float_512-256-180.raw";
const char writeFileName2[] = "sphere_singlepinhole_float_64-64-64.raw";

// void Projecion_SinglePinhole_3d(float* f,float* g);
// void backProjection_Singlepinhole_3d(float* f, float* g);
// void search_LessThan30_3d(float* f, float* g);
void knifeProjecion_SinglePinhole_3d(float* f, float* g, float* fov, Condition cond, int is_inverse = 0);
Eigen::Vector3f calculate_unit_vector(Eigen::Vector3f past, Eigen::Vector3f curr);



template <typename T>
void readRawFile (const char fname[], const size_t num, T* image);

template <typename T>
void writeRawFile (const char fname[], const size_t num, T* image);

int main(void)
{

	Condition cond;
	cond.img_w = 64;
	cond.img_h = 64;
	cond.img_d = 64;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 13;
	cond.thread_num = 4;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.2;
	cond.update_count = 1;
	cond.img_pixel_size = 0.2;
	cond.d_width = 0.08;
	cond.d_height = 0.08;

	float* f = (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));

	float* g = (float*)calloc(cond.detector_size_w * cond.detector_size_h * cond.detector_num, sizeof(float));
	float* fov = (float*)calloc(cond.detector_size_w * cond.detector_size_h, sizeof(float));
	float* h = (float*)calloc(cond.img_w * cond.img_h * cond.img_d, sizeof(float));

	readRawFile(readFileName, cond.img_w * cond.img_h * cond.img_d, f);
	readRawFile("fov_float_512-256.raw", cond.detector_size_w * cond.detector_size_h, fov);

	knifeProjecion_SinglePinhole_3d(f, g, fov, cond);
	writeRawFile(writeFileName1,  cond.detector_size_w * cond.detector_size_h * cond.detector_num, g);

	knifeProjecion_SinglePinhole_3d(g, h, fov, cond, 1);
	writeRawFile(writeFileName2,  cond.img_w * cond.img_h * cond.img_d, h);
}



void knifeProjecion_SinglePinhole_3d(float* f, float* g, float* fov, Condition cond, int is_inverse)
{
	int delta_detector = 360 / cond.detector_num;

	for(int theta_degree = 0; theta_degree < 360; theta_degree += delta_detector) {
		printf("theta_degree = %d\n", theta_degree);
		for(int m = 0; m < cond.detector_size_h; m++) {
			for(int n = 0; n < cond.detector_size_w; n++) {
				if(fov[m * cond.detector_size_w + n] == -1) { continue; }
				const float theta = theta_degree * M_PI / 180.0f;

				float r = cond.width_collimator / 2.;
				float collimator_x[7] = { 0., - sqrt(2. / 3.) * r, sqrt(2. / 3.) * r, - sqrt(1. / 6.) * r, - sqrt(1. / 6.) * r, sqrt(1. / 6.) * r, sqrt(1. / 6.) * r };
				float collimator_z[7] = { 0., 0., 0., - sqrt(1. / 6.) * r, sqrt(1. / 6.) * r, - sqrt(1. / 6.) * r, sqrt(1. / 6.) * r };

				for(int ray = 0; ray < 7; ray++)
				{
					float collimator_y = -1. * (cond.rotation_radius);

					Eigen::Vector3f on_collimator;
					on_collimator(0) = collimator_x[ray] * cosf(-theta) - collimator_y * sinf(-theta);
					on_collimator(1) = collimator_x[ray] * sinf(-theta) + collimator_y * cosf(-theta);
					on_collimator(2) = collimator_z[ray];

					float detector_x = (n - (cond.detector_size_w - 1.) / 2.) * cond.d_width;
					float detector_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);

					Eigen::Vector3f on_detector;
					on_detector(0) = detector_x * cosf(-theta) - detector_y * sinf(-theta);
					on_detector(1) = detector_x * sinf(-theta) + detector_y * cosf(-theta);
					on_detector(2) = ((cond.detector_size_h - 1.) / 2. - m) * cond.d_height;

					Eigen::Vector3f d = calculate_unit_vector(on_detector, on_collimator);

					//スケーリング処理
					// 多分on_collimatorは消してOK
					// on_collimator /= cond.img_pixel_size;
					Eigen::Vector3f sp = on_detector / cond.img_pixel_size;

					for (int sample_point = 0; sample_point < 300; sample_point++)
					{
						if(-(cond.img_w - 1.0f) / 2.0f < sp(0) && sp(0) < (cond.img_w - 1.0f) / 2.0f &&
						   -(cond.img_h - 1.0f) / 2.0f < sp(1) && sp(1) < (cond.img_h - 1.0f) / 2.0f &&
						   -(cond.img_d - 1.0f) / 2.0f < sp(2) && sp(2) < (cond.img_d - 1.0f) / 2.0f )
						{
							//双線形補完を行う左上の画素の座標(x0,y0)
							float x0 = floor(sp(0) - 0.5f) + 0.5f;
							float y0 = ceil(sp(1) - 0.5f) + 0.5f;
							float z0 = ceil(sp(2) - 0.5f) + 0.5f;

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

							float weight = (ray == 0) ? 1. / 4. : 1. / 8.;

							if(is_inverse)
							{
								int f_index = cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n;

								g[index1] += f[f_index] * V8 * weight;
								g[index2] += f[f_index] * V7 * weight;
								g[index3] += f[f_index] * V6 * weight;
								g[index4] += f[f_index] * V5 * weight;
								g[index5] += f[f_index] * V4 * weight;
								g[index6] += f[f_index] * V3 * weight;
								g[index7] += f[f_index] * V2 * weight;
								g[index8] += f[f_index] * V1 * weight;
							}
							else
							{
								float val = f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 + f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;
								val *= weight;

								g[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] += val;
							}
						}
						sp += d;
					}
				}
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


//
// void Projecion_SinglePinhole_3d(float* f,float* g)
// {
// 	//ピクセルサイズは0.17 cm × 0.17 cm
// 	float pixel_size = 0.17f;
// 	float* theta_collimator_xy = (float*)calloc( 512 * 256, sizeof(float));
// 	float* theta_collimator_zy = (float*)calloc( 512 * 256, sizeof(float));
//
// 	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy);
//
// 	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
// 	{
// 		cout << theta_degree << "° : projection processing..." << endl;
// 		const float theta = theta_degree * M_PI / 180.0f;
//
// 		//コリメータの座標の設定
// 		float collimator_x = 0.0f;
// 		float collimator_y = - 25.0f;
//
// 		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
// 		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
// 		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
// 		float z = 0.;
//
// 		for(int m = 0; m < 256; m++)
// 		{
// 			for(int n = 0; n < 512; n++)
// 			{
// 				//スケーリング処理
// 				float X = x / pixel_size;
// 				float Y = y / pixel_size;
// 				float Z = z / pixel_size;
//
// 				if(theta_collimator_xy[m * 512 + n] != 0.0f && theta_collimator_zy[m * 512 + n] != 0.0f)
// 				{
// 					for (int sp = 0; sp < 300; sp++)
// 					{
// 						if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f &&
// 						   -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f &&
// 						   -(L - 1.0f) / 2.0f < Z && Z < (L - 1.0f) / 2.0f )
// 						{
// 							//双線形補完を行う左上の画素の座標(x0,y0)
// 							float x0 = floor(X - 0.5f) + 0.5f;
// 							float y0 = ceil(Y - 0.5f) + 0.5f;
// 							float z0 = ceil(Z - 0.5f) + 0.5f;
//
// 							//i,jに戻す
// 							float J = x0 + (W - 1.0f) / 2.0f;
// 							float I = (H - 1.0f) / 2.0f - y0;
// 							float K = (L - 1.0f) / 2.0f - z0;
//
// 							int index1 = W * H * K + W * I + J;
// 							int index2 = W * H * K + W * I + (J + 1);
// 							int index3 = W * H * K + W * (I + 1) + J;
// 							int index4 = W * H * K + W * (I + 1) + (J + 1);
// 							int index5 = W * H * (K + 1) + W * I + J;
// 							int index6 = W * H * (K + 1) + W * I + (J + 1);
// 							int index7 = W * H * (K + 1) + W * (I + 1) + J;
// 							int index8 = W * H * (K + 1) + W * (I + 1) + (J + 1);
//
// 							//(X , Y)から左上までの距離dx, dy, dz
// 							float dx = fabs(X - x0);
// 							float dy = fabs(Y - y0);
// 							float dz = fabs(Z - z0);
//
// 							float V1 = dx * dy * dz;
// 							float V2 = (1.0f - dx) * dy * dz;
// 							float V3 = dx * (1.0f - dy) * dz;
// 							float V4 = (1.0f - dx) * (1.0f - dy) * dz;
// 							float V5 = dx * dy * (1.0f - dz);
// 							float V6 = (1.0f - dx) * dy * (1.0f - dz);
// 							float V7 = dx * (1.0f - dy) * (1.0f - dz);
// 							float V8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);
//
// 							g[512 * 256 * theta_degree + 512 * m + n] += f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 +
// 							 									         f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;
//
// 						}
//
// 						X +=  cos(theta_collimator_zy[512 * m + n]) * sin(theta - theta_collimator_xy[512 * m + n]);
// 						Y +=  cos(theta_collimator_zy[512 * m + n]) * cos(theta - theta_collimator_xy[512 * m + n]);
// 						Z +=  -sin(theta_collimator_zy[512 * m + n]);
// 					}
// 				}
// 			}
// 		}
// 	}
// }
//
//
//
// void backProjection_Singlepinhole_3d(float* f, float* g)
// {
// 	//ピクセルサイズは0.17 cm × 0.17 cm
// 	float pixel_size = 0.17f;
// 	float* theta_collimator_xy = (float*)calloc( 512 * 256, sizeof(float));
// 	float* theta_collimator_zy = (float*)calloc( 512 * 256, sizeof(float));
//
// 	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy);
//
// 	for(int theta_degree = 0; theta_degree < 360; theta_degree++)
// 	{
// 		cout << theta_degree << "° : backprojection processing..." << endl;
// 		const float theta = theta_degree * M_PI / 180.0f;
//
// 		//コリメータの座標の設定
// 		float collimator_x = 0.0f;
// 		float collimator_y = -25.0f;
//
// 		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
// 		float x = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
// 		float y = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
// 		float z = 0.;
//
// 		for(int m = 0; m < 256; m++)
// 		{
// 			for(int n = 0; n < 512; n++)
// 			{
// 				//スケーリング処理
// 				float X = x / pixel_size;
// 				float Y = y / pixel_size;
// 				float Z = z / pixel_size;
//
//
// 				if(theta_collimator_xy[m * 512 + n] != 0.0f && theta_collimator_zy[m * 512 + n] != 0.0f)
// 				{
// 					for (int sp = 0; sp < 300; sp++)
// 					{
// 						if(-(W - 1.0f) / 2.0f < X && X < (W - 1.0f) / 2.0f &&
// 						   -(H - 1.0f) / 2.0f < Y && Y < (H - 1.0f) / 2.0f &&
// 						   -(L - 1.0f) / 2.0f < Z && Z < (L - 1.0f) / 2.0f )
// 						{
// 							//双線形補完を行う左上の画素の座標(x0,y0)
// 							float x0 = floor(X - 0.5f) + 0.5f;
// 							float y0 = ceil(Y - 0.5f) + 0.5f;
// 							float z0 = ceil(Z - 0.5f) + 0.5f;
//
// 							//i,jに戻す
// 							float J = x0 + (W - 1.0f) / 2.0f;
// 							float I = (H - 1.0f) / 2.0f - y0;
// 							float K = (L - 1.0f) / 2.0f - z0;
//
// 							int index1 = W * H * K + W * I + J;
// 							int index2 = W * H * K + W * I + (J + 1);
// 							int index3 = W * H * K + W * (I + 1) + J;
// 							int index4 = W * H * K + W * (I + 1) + (J + 1);
// 							int index5 = W * H * (K + 1) + W * I + J;
// 							int index6 = W * H * (K + 1) + W * I + (J + 1);
// 							int index7 = W * H * (K + 1) + W * (I + 1) + J;
// 							int index8 = W * H * (K + 1) + W * (I + 1) + (J + 1);
//
// 							//(X , Y)から左上までの距離dx, dy, dz
// 							float dx = fabs(X - x0);
// 							float dy = fabs(Y - y0);
// 							float dz = fabs(Z - z0);
//
// 							float V1 = dx * dy * dz;
// 							float V2 = (1.0f - dx) * dy * dz;
// 							float V3 = dx * (1.0f - dy) * dz;
// 							float V4 = (1.0f - dx) * (1.0f - dy) * dz;
// 							float V5 = dx * dy * (1.0f - dz);
// 							float V6 = (1.0f - dx) * dy * (1.0f - dz);
// 							float V7 = dx * (1.0f - dy) * (1.0f - dz);
// 							float V8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);
//
//
// 							g[index1] += f[512 * 256 * theta_degree + 512 * m + n] * V8;
// 							g[index2] += f[512 * 256 * theta_degree + 512 * m + n] * V7;
// 							g[index3] += f[512 * 256 * theta_degree + 512 * m + n] * V6;
// 							g[index4] += f[512 * 256 * theta_degree + 512 * m + n] * V5;
// 							g[index5] += f[512 * 256 * theta_degree + 512 * m + n] * V4;
// 							g[index6] += f[512 * 256 * theta_degree + 512 * m + n] * V3;
// 							g[index7] += f[512 * 256 * theta_degree + 512 * m + n] * V2;
// 							g[index8] += f[512 * 256 * theta_degree + 512 * m + n] * V1;
// 						}
//
// 						X +=  cos(theta_collimator_zy[512 * m + n]) * sin(theta - theta_collimator_xy[512 * m + n]);
// 						Y +=  cos(theta_collimator_zy[512 * m + n]) * cos(theta - theta_collimator_xy[512 * m + n]);
// 						Z +=  -sin(theta_collimator_zy[512 * m + n]);
// 					}
// 				}
// 			}
// 		}
// 	}
// }
//
// void search_LessThan30_3d(float* f, float* g)
// {
// 	//確認用
// 	float* h = (float*)calloc(512 * 256, sizeof(float));
//
// 	for(int i = 0; i < 256; i++)
// 	{
// 		for(int j = 0; j < 512; j++)
// 		{
// 			float i0 = i * 0.08f + 0.04f;
// 			float j0 = j * 0.08f + 0.04f;
// 			float x = - 512 * 0.08f / 2 + j0;
// 			float z =   256 * 0.08f / 2 - i0;
// 			float theta_Ditector_xy = atanf(x / 7.5f);
// 			float theta_Ditector_zy = atanf(z / 7.5f);
//
// 			// 半径7.5/sqrt(3)はピンホールに対する角度が30の時の円の半径
// 			float radius = 7.5 / sqrt(3);
//
// 			if(pow(x, 2.0f) + pow(z, 2.0f) < pow(radius, 2.0f))
// 			{
// 				f[i * 512 + j] = theta_Ditector_xy;
// 				g[i * 512 + j] = theta_Ditector_zy;
// 				h[i * 512 + j] = 100.0f;//確認用
// 			}
// 			else
// 			{
// 				f[i * 512 + j] = 0.0f;
// 				g[i * 512 + j] = 0.0f;
// 				h[i * 512 + j] = 0.;//確認用
// 			}
// 		}
// 	}
// 	writeRawFile("circleLessThan30_float_512-256.raw", 512 * 256, h);
// }


template <typename T>
void readRawFile (const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
}


template <typename T>
void writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}
