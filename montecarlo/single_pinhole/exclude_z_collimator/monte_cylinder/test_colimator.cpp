
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <random>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define detector_size 65.
#define mumap_size 128

int size = (int)detector_size;


int main()
{
  // 条件
  float rotation_radius = 15.;
  float distance_collimator_to_detector = 7.5;
  float height_collimator = 1.;
  float width_collimator = 0.3;

  Eigen::Vector2f past_rotated;
  Eigen::Vector2f curr_rotated;

  past_rotated(0) = 14.5;
  past_rotated(1) = (0.28867513 + 0.15);
  curr_rotated(0) = 15.5;
  curr_rotated(1) = -(0.28867513 + 0.15);

  if(abs(curr_rotated(0) - past_rotated(0)) < 0.00001) { exit(-1); }

  Eigen::Vector2f photon_vec;
  photon_vec << curr_rotated(0) - past_rotated(0), curr_rotated(1) - past_rotated(1);

  if(abs(photon_vec(0)) < 0.00001) { exit(-1); }
  float tan_collimator = tan(photon_vec(1) / photon_vec(0));

  // コリメータ手前
  float vec_scale = (rotation_radius - height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
  float y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
  if(abs(y_on_collimator) > width_collimator / 2. + height_collimator / 2. * tan(M_PI / 6)) { exit(-1); }

  // コリメータ奥
  vec_scale = (rotation_radius + height_collimator / 2 - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
  y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
  if(abs(y_on_collimator) > width_collimator / 2. + height_collimator / 2. * tan(M_PI / 6)) { exit(-1); }

  // コリメータ真ん中
  vec_scale = (rotation_radius - past_rotated(0)) / (curr_rotated(0) - past_rotated(0));
  y_on_collimator = past_rotated(1) + vec_scale * (curr_rotated(1) - past_rotated(1));
  if(abs(y_on_collimator) > width_collimator / 2.) { exit(-1); }


  //検出器のピクセルサイズが0.5 cmであるため ×2をしている
  // yの値が大きい→検出器の番号は小さい
  float y_on_detector = y_on_collimator - distance_collimator_to_detector * tan_collimator * 2;

  float i0 = size / 2. + y_on_detector;

  if(i0 < 0. || i0 > 65.) { exit(-1); }

  int i1 = (int)floor(i0);

  cout << "tan_collimator = " << tan_collimator << endl;
  cout << "y_on_detector = " << y_on_detector << endl;
  cout << "i0 = " << i0 << endl;
  cout << "i1 = " << i1 << endl;
}
