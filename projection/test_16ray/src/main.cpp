/*

感度補正~ 再構成まで

floatにしてある
genrand_real1() //一様実乱数[0,1] (32ビット精度)
genrand_real2() //一様実乱数[0,1) (32ビット精度)
genrand_real3() //一様実乱数(0,1) (32ビット精度)

cond.thread_num = 1にすると再構成うまくいく

FIXME : 7rayの半径のとり方が円と同じになってる

FIXME : 今primaryしか投影データとってない
*/

#include <stdio.h>
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
// #define DEBUG


enum PinholeShape { CIRCLE, RECTANGLE };

typedef struct {
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int surface_source_w;
	int surface_source_h;
	int img_w;
	int img_h;
	int img_d;
	int update_count;
	int photon_scale;
	int photon_num_efficiency_map;
	int max_scattering_count;
	int pinhole_count;
	int pinhole_img_w;
	int pinhole_img_h;
	int pinhole_shape;
	int thread_num;
	int ray_num;
	bool is_voxel;
	bool has_efficiency_map;
	float cut_off_energy;
	float rotation_radius;
	float distance_collimator_to_detector;
	float collimator_h;
	/*------------------*/
	// 円形ピンホールの場合
	float collimator_w;
	float aperture_degree;
	// 矩形ピンホールの場合
	float collimator_w_y_axis;
	float collimator_w_z_axis;
	float aperture_degree_xy;
	float aperture_degree_zx;
	/*------------------*/
	float img_pixel_size;
	float detector_pixel_size_w;
	float detector_pixel_size_h;
	float pinhole_img_pixel_size;
	float time;
} Condition;


void launch_test_16ray(Condition cond);
void projecion_multi_pinhole_3d(vector<float>& f, vector<float>& g, vector<float>& absorp_map, vector<int> fov, vector<float>& pinhole_theta_xy, vector<float>& pinhole_theta_zx, vector<float>& pinhole_center, vector<float>& collimator_x, vector<float>& collimator_z, vector<float>& pass_collimator, Condition cond, int is_inverse = 0);
void set_projection_line_on_collimator(vector<float> pass_collimator, vector<float> pinhole_theta_xy, vector<float> pinhole_theta_zx, vector<float> pinhole_center, vector<float> x, vector<float> z, int pinhole_num, Condition cond);
void rotate_axis(Eigen::Vector3f &v, float theta);


int main()
{
	Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.img_d = 128;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.surface_source_w = 2048;
	cond.surface_source_h = 1024;
	cond.rotation_radius = 25;
	// cond.photon_scale = 15000;
	cond.photon_scale = 15000;
	cond.photon_num_efficiency_map = 1E6;
	cond.pinhole_img_w = 2048;
	cond.pinhole_img_h = 1024;
	cond.thread_num = 4;
	cond.ray_num = 16;
	cond.max_scattering_count = 5;
	cond.cut_off_energy = 30.;
	cond.distance_collimator_to_detector = 7.6;
	cond.collimator_h = 1.;
	/*-----------------------------*/
	// 円形ピンホールの場合
	cond.collimator_w = 0.5;
	cond.aperture_degree = 24.0f;
	// 矩形ピンホールの場合
	cond.collimator_w_y_axis = 0.5;
	cond.collimator_w_z_axis = 0.5;
	cond.aperture_degree_xy = 24.;
	cond.aperture_degree_zx = 24.;
	/*-----------------------------*/
	cond.pinhole_shape = RECTANGLE;
	cond.update_count = 100;
	cond.img_pixel_size = 0.16;
	cond.detector_pixel_size_w = 0.08;
	cond.detector_pixel_size_h = 0.08;
	cond.pinhole_img_pixel_size = 0.02;
	cond.has_efficiency_map = true;

	

	launch_test_16ray(cond);
}

void launch_test_16ray(Condition cond)
{
	// 全部1の配列用意

	// これを16rayで投影
	// 書き込み
	vector<float> f(cond.img_w * cond.img_h * cond.img_d, 1.0f);
	vector<float> f2(cond.img_w * cond.img_h * cond.img_d, 0.0f);
	vector<float> g(cond.detector_size_w * cond.detector_size_h * cond.detector_num, 0.0f);
	float r_x = cond.collimator_w / 2.0f;
	float x1 = - r_x * 3. / 4;
	float x2 =  - r_x * 1. / 4;
	float x3 =  r_x * 1. / 4;
	float x4 =  r_x * 3. / 4;
	float r_y = cond.collimator_h / 2.0f;
	float y1 = r_y * 3. / 4;
	float y2 = r_y * 1. / 4;
	float y3 = - r_y * 1. / 4;
	float y4 = - r_y * 3. / 4;

	std::vector<float> pinhole_theta_xy{ -21.56, 21.56, 0., -21.56, 21.56 };
	cond.pinhole_count = pinhole_theta_xy.size();
	std::vector<float> pinhole_theta_zx{ 	-8.35f, -8.35f, 0., 8.35f, 8.35f };
	std::vector<float> pinhole_center
	{
	    cond.rotation_radius, - 9.,  3.62,
	    cond.rotation_radius,   9.,  3.62,
	    cond.rotation_radius,    0.,  0.,
	    cond.rotation_radius, - 9., - 3.62,
	    cond.rotation_radius,   9., - 3.62
	};
	std::vector<float> collimator_x { x1, x2, x3, x4, x1, x2, x3, x4, x1, x2, x3, x4, x1, x2, x3, x4 };
	std::vector<float> collimator_z { y1, y2, y3, y4, y1, y2, y3, y4, y1, y2, y3, y4, y1, y2, y3, y4 };
	std::vector<float> pass_collimator(coord_num * cond.ray_num, 0.);
	vector<int> fov(cond.detector_size_w * cond.detector_size_h);
	readRawFile("img/fov_5pinhole_square_5-5mm_int_512-256.raw", fov);
	projecion_multi_pinhole_3d(f, g, f2, fov, pinhole_theta_xy, pinhole_theta_zx, pinhole_center, collimator_x, collimator_z, pass_collimator, cond);
	writeRawFile("result/f_proj_float_512-256-180.raw", g);
}

// FIXME : 7raysにすると画像の下の方がおかしくなる
void projecion_multi_pinhole_3d(vector<float>& f, vector<float>& g, vector<float>& absorp_map, vector<int> fov, vector<float>& pinhole_theta_xy, vector<float>& pinhole_theta_zx, vector<float>& pinhole_center, vector<float>& collimator_x, vector<float>& collimator_z, vector<float>& pass_collimator, Condition cond, int is_inverse)
{
	int delta_detector = 360 / cond.detector_num;
	// rep(idx, cond.detector_num)
	rep(idx, 1)
	{
		int theta_degree = (idx) * delta_detector;
		rep(m, cond.detector_size_h)rep(n, cond.detector_size_w)
		{
			int pinhole_num = fov[m * cond.detector_size_w + n];
			if(pinhole_num == -1) continue;
			const float theta = theta_degree * M_PI / 180.0f;

			set_projection_line_on_collimator(pass_collimator, pinhole_theta_xy, pinhole_theta_zx, pinhole_center, collimator_x, collimator_z, pinhole_num, cond);

			for(int ray = 0; ray < cond.ray_num; ray++)
			{
				int count_h2o = 0;
				int count_ca = 0;

				Eigen::Vector3f on_collimator;
				on_collimator(0) = pass_collimator[ray * coord_num + X] * cosf(-theta) - pass_collimator[ray * coord_num + Y] * sinf(-theta);
				on_collimator(1) = pass_collimator[ray * coord_num + X] * sinf(-theta) + pass_collimator[ray * coord_num + Y] * cosf(-theta);
				on_collimator(2) = pass_collimator[ray * coord_num + Z];


				// if(theta_degree == 0 && ray == 0 && n == cond.detector_size_w / 2 && m == cond.detector_size_h / 2 && is_inverse == 0)
				// {
				// 	printf("on_collimator(0) %f\n", on_collimator(0));
				// 	printf("on_collimator(1) %f\n", on_collimator(1));
				// 	printf("on_collimator(2) %f\n", on_collimator(2));
				// }

				float detector_x = (n - (cond.detector_size_w - 1.) / 2.) * cond.detector_pixel_size_w;
				float detector_y = -1. * (cond.rotation_radius + cond.distance_collimator_to_detector);

				Eigen::Vector3f on_detector;
				on_detector(0) = detector_x * cosf(-theta) - detector_y * sinf(-theta);
				on_detector(1) = detector_x * sinf(-theta) + detector_y * cosf(-theta);
				on_detector(2) = ((cond.detector_size_h - 1.) / 2. - m) * cond.detector_pixel_size_h;

				Eigen::Vector3f direction_sp = calculate_unit_vector(on_detector, on_collimator);
				// printf("direction_sp(X) = %f\n", direction_sp(X));
				// printf("direction_sp(Y) = %f\n", direction_sp(Y));
				// printf("direction_sp(Z) = %f\n", direction_sp(Z));

				//スケーリング処理
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

						// if(theta_degree == 0 && ray == 0 && n == cond.detector_size_w / 2 && m == cond.detector_size_h / 2 && sample_point > 70 && sample_point < 140 && is_inverse == 0)
						// {
						// 	printf("sample_point =  %d\n", sample_point);
						// 	printf("sp(0) =  %f\n", sp(0));
						// 	printf("sp(1) =  %f\n", sp(1));
						// 	printf("sp(2) =  %f\n", sp(2));
						// 	printf("x0 =  %f\n", x0);
						// 	printf("x0 =  %f\n", x0);
						// 	printf("y0 =  %f\n", y0);
						// 	printf("z0 =  %f\n", z0);
						// 	printf("J =  %f\n", J);
						// 	printf("I =  %f\n", I);
						// 	printf("D =  %f\n", D);
						// }

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

						if(is_inverse)
						{
							int f_index = cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n;
							
							g[index1] += f[f_index] * V8;
							g[index2] += f[f_index] * V7;
							g[index3] += f[f_index] * V6;
							g[index4] += f[f_index] * V5;
							g[index5] += f[f_index] * V4;
							g[index6] += f[f_index] * V3;
							g[index7] += f[f_index] * V2;
							g[index8] += f[f_index] * V1;
						}
						else
						{
							float val = f[index1] * V8 + f[index2] * V7 + f[index3] * V6 + f[index4] * V5 + f[index5] * V4 + f[index6] * V3 + f[index7] * V2 + f[index8] * V1;

							/* ここってimg_pixel_sizeかけなきゃいけないの？？ */
							// int j = round(sp(0) * cond.img_pixel_size + (cond.img_w - 1.0f) / 2.0f);
							// int i = round((cond.img_h - 1.0f) / 2.0f - sp(1) * cond.img_pixel_size);
							// int d = round((cond.img_d - 1.0f) / 2.0f - sp(2) * cond.img_pixel_size);
							// if(absorp_map[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.26) { count_ca++; }
							// if(absorp_map[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.14) { count_h2o++; }

							// if(theta_degree == 0 && ray == 0 && n == cond.detector_size_w / 2 && m == cond.detector_size_h / 2 && sample_point > 70 && sample_point < 140 && is_inverse == 0)
							// {
							// 	printf("i = %d\n", i);
							// 	printf("j = %d\n", j);
							// 	printf("d = %d\n", d);
							//
							// 	int j1 = sp(0) + (cond.img_w - 1) / 2;
							// 	int i1 = (cond.img_h - 1) / 2 - sp(1);
							// 	int d1 = (cond.img_d - 1) / 2 - sp(2);
							//
							// 	printf("i1 = %d\n", i1);
							// 	printf("j1 = %d\n", j1);
							// 	printf("d1 = %d\n", d1);
							// }

							int j = round(sp(0) + (cond.img_w - 1.) / 2.);
							int i = round((cond.img_h - 1.) / 2. - sp(1));
							int d = round((cond.img_d - 1.) / 2. - sp(2));
							if(absorp_map[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.26f) { count_ca++; }
							if(absorp_map[d * cond.img_w * cond.img_h + i * cond.img_w + j] > 0.14f) { count_h2o++; }

							// if(theta_degree == 0 && ray == 0 && n == cond.detector_size_w / 2 && m == cond.detector_size_h / 2 && sample_point > 70 && sample_point < 140 && is_inverse == 0)
							// {
							// 	printf("count_ca = %d\n", count_ca);
							// 		printf("count_h2o = %d\n", count_h2o);
							// }

							float mu_h2o = 0.1538;
							float mu_ca = 0.2748;
							g[cond.detector_size_w * cond.detector_size_h * (theta_degree / delta_detector) + cond.detector_size_w * m + n] += val * exp(- mu_h2o * count_h2o * cond.img_pixel_size - mu_ca * count_ca * cond.img_pixel_size);
						}
					}
					sp += direction_sp;
				}
			}
		}
	}
}


void rotate_axis(Eigen::Vector3f &v, float theta)
{
	Eigen::Vector3f rotated;
	rotated(0) = v(0) * cos(theta) - v(1) * sin(theta);
	rotated(1) = v(0) * sin(theta) + v(1) * cos(theta);
	rotated(2) = v(2);

	v = rotated;
}

void set_projection_line_on_collimator(vector<float> pass_collimator, vector<float> pinhole_theta_xy, vector<float> pinhole_theta_zx, vector<float> pinhole_center, vector<float> x, vector<float> z, int pinhole_num, Condition cond)
{
	// float y = -1. * (cond.rotation_radius);

	Eigen::Vector3f pinhole_ctr; pinhole_ctr << pinhole_center[pinhole_num * coord_num + X], pinhole_center[pinhole_num * coord_num + Y], pinhole_center[pinhole_num * coord_num + Z];
	float theta_xy = pinhole_theta_xy[pinhole_num];
	float theta_zx = pinhole_theta_zx[pinhole_num];
	for(int i = 0; i < cond.ray_num; i++)
	{
		float degree_to_rad = M_PI / 180.;
		Eigen::Vector3f on_collimator;  on_collimator << x[i], 0., z[i];

		// 90度回転（検出器の位置をx軸正にする）
		float theta_90 = M_PI / 2.;
		rotate_axis(on_collimator, theta_90);

		Eigen::Vector3f on_collimator_rot;
		{
			Eigen::Matrix3f rot_xy, rot_zx;
			Eigen::Vector3f axis; axis << 0, 0, 1;
			rot_xy = Eigen::AngleAxisf( theta_xy * degree_to_rad, axis );

			axis << 0, 1, 0;
			rot_zx = Eigen::AngleAxisf( theta_zx * degree_to_rad, axis );

			on_collimator_rot = rot_zx * rot_xy * on_collimator + pinhole_ctr;
		}
		// pass_collimator[i] << on_collimator_rot(X), on_collimator_rot(Y), on_collimator_rot(Z);

		// -90度回転（検出器の位置をy軸負に戻す）
		rotate_axis(on_collimator_rot, -theta_90);

		pass_collimator[i * coord_num + X] = on_collimator_rot(X);
		pass_collimator[i * coord_num + Y] = on_collimator_rot(Y);
		pass_collimator[i * coord_num + Z] = on_collimator_rot(Z);
	}
}
