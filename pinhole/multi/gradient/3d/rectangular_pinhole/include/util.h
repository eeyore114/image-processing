#ifndef  _UTIL_H_
#define  _UTIL_H_

#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/Geometry"

bool is_out_of_image(int i, int height, int j, int width);
bool is_in_image(int i, int height, int j, int width);
bool is_out_of_image_3d(int i, int height, int j, int width, int k, int depth);
bool is_in_image_3d(int i, int height, int j, int width, int k, int depth);
Eigen::Vector3f calculate_unit_vector(Eigen::Vector3f past, Eigen::Vector3f curr);
Eigen::Vector2f calculate_unit_vector_2f(Eigen::Vector2f past, Eigen::Vector2f curr);
int transform_img_coordinate_same_axis(float x, int w, float pixel_size = 1.);
int transform_img_coordinate_opposite_axis(float y, int h, float pixel_size = 1.);
float transform_origin_coordinate_same_axis(int j, int w, float pixel_size = 1.);
float transform_origin_coordinate_opposite_axis(int i, int h, float pixel_size = 1.);
float degree_to_rad(float theta_degree);
float rad_to_degree(float theta);



#endif /*_UTIL_H_*/
