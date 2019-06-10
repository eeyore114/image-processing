/*
4/29 19:19
search_LessThan30_3d作成中

コリメータ通るかどうかの判断の部分作成中

4/30

とりあえず作成してみた。まだ確認はしていない。
コリメータの座標を決めている場所がよくわかんない


一応できてそう。voxelから1000個の投影画像で再構成したら円の形は出てきた。
1万か10万の光子数でやったものを確かめてみる。
*/

#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"


template <class T>
void readRawFile (string fname, const size_t num, T* image);

template <class T>
void writeRawFile (string fname, const size_t num, T* image);

void backProjection_Singlepinhole_3d();
void search_LessThan30_3d(float* f, float* g, int detector_size_w, int detector_size_h, float distance_collimator_to_detector, float rotation_radius, float height_collimator, float width_collimator);


int main()
{
	backProjection_Singlepinhole_3d();
}


void backProjection_Singlepinhole_3d()
{
	string read_file_name = "detector_single_pinhole_float_180-180-180.raw";
	// string read_file_name = "Brain_3d_singlepinhole_float_512-256-360.raw";
	string write_file_name = "backproj_float_128-128-128.raw";
	float rotation_radius = 10.;
	float distance_collimator_to_detector = 7.5;
	float height_collimator = 1.;
	float width_collimator = 0.5;

	const int MUMAP_SIZE_W = 128;
	const int MUMAP_SIZE_H = 128;
	const int MUMAP_SIZE_D = 128;
	const int DETECTOR_NUM = 180;
	const int DETECTOR_SIZE_W = 180;
	const int DETECTOR_SIZE_H = 180;

	float* f = (float*)calloc(DETECTOR_NUM * DETECTOR_SIZE_W * DETECTOR_SIZE_H, sizeof(float));
	float* g = (float*)calloc(MUMAP_SIZE_W * MUMAP_SIZE_H * MUMAP_SIZE_D, sizeof(float));
	readRawFile(read_file_name, DETECTOR_NUM * DETECTOR_SIZE_W * DETECTOR_SIZE_H, f);

	//ピクセルサイズは0.17 cm × 0.17 cm
	// float pixel_size = 0.17f;

	float pixel_size = 0.35;
	float* theta_collimator_xy = (float*)calloc(DETECTOR_SIZE_W * DETECTOR_SIZE_H, sizeof(float));
	float* theta_collimator_zy = (float*)calloc(DETECTOR_SIZE_W * DETECTOR_SIZE_H, sizeof(float));

	search_LessThan30_3d(theta_collimator_xy, theta_collimator_zy, DETECTOR_SIZE_W, DETECTOR_SIZE_H, distance_collimator_to_detector, rotation_radius, height_collimator, width_collimator);

	for(int theta_degree = 0; theta_degree < 360; theta_degree += 2)
	{
		cout << theta_degree << "° : backprojection processing..." << endl;
		const float theta = theta_degree * M_PI / 180.0f;

		//コリメータの座標の設定
		float collimator_x = 0.0f;
		float collimator_y = -1. * (rotation_radius + distance_collimator_to_detector);


		//コリメータの座標を原点を中心に-theta回転させた座標(最初のサンプルポイント)
		Eigen::Vector3f on_collimator;
		on_collimator(0) = collimator_x * cosf(-theta) - collimator_y * sinf(-theta);
		on_collimator(1) = collimator_x * sinf(-theta) + collimator_y * cosf(-theta);
		on_collimator(2) = 0.;

		for(int m = 0; m < DETECTOR_SIZE_H; m++)
		{
			for(int n = 0; n < DETECTOR_SIZE_W; n++)
			{
				//スケーリング処理
				float X = on_collimator(0) / pixel_size;
				float Y = on_collimator(1) / pixel_size;
				float Z = on_collimator(2) / pixel_size;


				if(theta_collimator_xy[m * DETECTOR_SIZE_W + n] != 0.0f && theta_collimator_zy[m * DETECTOR_SIZE_W + n] != 0.0f)
				{
					for (int sp = 0; sp < 300; sp++)
					{
						if(-(MUMAP_SIZE_W - 1.0f) / 2.0f < X && X < (MUMAP_SIZE_W - 1.0f) / 2.0f &&
						   -(MUMAP_SIZE_H - 1.0f) / 2.0f < Y && Y < (MUMAP_SIZE_H - 1.0f) / 2.0f &&
						   -(MUMAP_SIZE_D - 1.0f) / 2.0f < Z && Z < (MUMAP_SIZE_D - 1.0f) / 2.0f )
						{
							//双線形補完を行う左上の画素の座標(x0,y0)
							float x0 = floor(X - 0.5f) + 0.5f;
							float y0 = ceil(Y - 0.5f) + 0.5f;
							float z0 = ceil(Z - 0.5f) + 0.5f;

							//i,jに戻す
							float J = x0 + (MUMAP_SIZE_W - 1.0f) / 2.0f;
							float I = (MUMAP_SIZE_H - 1.0f) / 2.0f - y0;
							float D = (MUMAP_SIZE_D - 1.0f) / 2.0f - z0;

							int index1 = MUMAP_SIZE_W * MUMAP_SIZE_H * D + MUMAP_SIZE_W * I + J;
							int index2 = MUMAP_SIZE_W * MUMAP_SIZE_H * D + MUMAP_SIZE_W * I + (J + 1);
							int index3 = MUMAP_SIZE_W * MUMAP_SIZE_H * D + MUMAP_SIZE_W * (I + 1) + J;
							int index4 = MUMAP_SIZE_W * MUMAP_SIZE_H * D + MUMAP_SIZE_W * (I + 1) + (J + 1);
							int index5 = MUMAP_SIZE_W * MUMAP_SIZE_H * (D + 1) + MUMAP_SIZE_W * I + J;
							int index6 = MUMAP_SIZE_W * MUMAP_SIZE_H * (D + 1) + MUMAP_SIZE_W * I + (J + 1);
							int index7 = MUMAP_SIZE_W * MUMAP_SIZE_H * (D + 1) + MUMAP_SIZE_W * (I + 1) + J;
							int index8 = MUMAP_SIZE_W * MUMAP_SIZE_H * (D + 1) + MUMAP_SIZE_W * (I + 1) + (J + 1);

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

							g[index1] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V8;
							g[index2] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V7;
							g[index3] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V6;
							g[index4] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V5;
							g[index5] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V4;
							g[index6] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V3;
							g[index7] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V2;
							g[index8] += f[DETECTOR_SIZE_W * DETECTOR_SIZE_H * (theta_degree / 2) + DETECTOR_SIZE_W * m + n] * V1;
						}

						X +=  cos(theta_collimator_zy[DETECTOR_SIZE_W * m + n]) * sin(theta - theta_collimator_xy[DETECTOR_SIZE_W * m + n]);
						Y +=  cos(theta_collimator_zy[DETECTOR_SIZE_W * m + n]) * cos(theta - theta_collimator_xy[DETECTOR_SIZE_W * m + n]);
						Z +=  -sin(theta_collimator_zy[DETECTOR_SIZE_W * m + n]);
					}
				}
			}
		}
	}


	writeRawFile(write_file_name, MUMAP_SIZE_W * MUMAP_SIZE_H * MUMAP_SIZE_D, g);
}


void search_LessThan30_3d(float* f, float* g, int detector_size_w, int detector_size_h, float distance_collimator_to_detector, float rotation_radius, float height_collimator, float width_collimator)
{
	//確認用
	float* h = (float*)calloc(detector_size_w * detector_size_h, sizeof(float));

	// 検出器の幅
	float d_width = 0.5;
	float d_height = 0.5;
	for(int i = 0; i < detector_size_h; i++)
	{
		for(int j = 0; j < detector_size_w; j++)
		{
			float i0 = i * d_width + d_width / 2.;
			float j0 = j * d_height + d_height / 2.;
			Eigen::Vector3f on_detector;
			on_detector(0) = - detector_size_w * d_width / 2. + j0;
			on_detector(1) = -1. * (rotation_radius + distance_collimator_to_detector);
			on_detector(2) =   detector_size_h * d_height / 2. - i0;
			float theta_Ditector_xy = atanf(on_detector(0) / distance_collimator_to_detector);
			float theta_Ditector_zy = atanf(on_detector(2) / distance_collimator_to_detector);

			Eigen::Vector3f center_colimator;
			center_colimator << 0., -1. * rotation_radius, 0.;

			// コリメータ奥
			Eigen::Vector3f on_colimator;
			on_colimator(1) = -1. * (rotation_radius + height_collimator / 2.);
			float vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
			on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
			on_colimator(2) = on_detector(2) + vec_scale * (center_colimator(2) - on_detector(2));
			float colimator_radius = width_collimator / 2. + tan(M_PI / 6);

			if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.))
			{
				f[i * detector_size_w + j] = 0.0f;
				g[i * detector_size_w + j] = 0.0f;
				h[i * detector_size_w + j] = 0.;//確認用
				continue;
			}

			// コリメータ手前
			on_colimator(1) = -1. * (rotation_radius - height_collimator / 2.);
			vec_scale = (on_colimator(1) - on_detector(1)) / (center_colimator(1) - on_detector(1));
			on_colimator(0) = on_detector(0) + vec_scale * (center_colimator(0) - on_detector(0));
			on_colimator(2) = on_detector(2) + vec_scale * (center_colimator(2) - on_detector(2));
			colimator_radius = width_collimator / 2. + tan(M_PI / 6);

			if(pow(on_colimator(0), 2.) + pow(on_colimator(2), 2.) > pow(colimator_radius, 2.))
			{
				f[i * detector_size_w + j] = 0.0f;
				g[i * detector_size_w + j] = 0.0f;
				h[i * detector_size_w + j] = 0.;//確認用
				continue;
			}

			f[i * detector_size_w + j] = theta_Ditector_xy;
			g[i * detector_size_w + j] = theta_Ditector_zy;
			h[i * detector_size_w + j] = 100.0f;//確認用
		}
	}
	string write_file_name = "test_within_30_float_180-180.raw";
	writeRawFile(write_file_name, detector_size_w * detector_size_h, h);
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
