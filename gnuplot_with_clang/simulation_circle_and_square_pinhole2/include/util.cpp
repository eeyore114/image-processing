#include "util.h"
#define rep(i,n) for(int (i)=0;(i)<(n);(i)++)

enum Coordinate { X, Y, Z, coord_num };

bool is_in_image(int i, int height, int j, int width) { return (0 < i && i < height && 0 < j && j < width); }
bool is_out_of_image(int i, int height, int j, int width) { return (0 > i || i > height || 0 > j || j > width); }
bool is_in_image_3d(int i, int height, int j, int width, int k, int depth) { return (0 < i && i < height && 0 < j && j < width && 0 < k && k < depth); }
bool is_out_of_image_3d(int i, int height, int j, int width, int k, int depth) { return (0 > i || i > height || 0 > j || j > width || 0 > k || k > depth); }
int transform_img_coordinate_same_axis(float x, int w, float pixel_size) { return x / pixel_size + (w - 1.0f) / 2.0f; }
int transform_img_coordinate_opposite_axis(float y, int h, float pixel_size) { return (h - 1.0f) / 2.0f - y / pixel_size; }
float transform_origin_coordinate_same_axis(int j, int w, float pixel_size) { return  (j - (w - 1.0f) / 2.0f) * pixel_size; }
float transform_origin_coordinate_opposite_axis(int i, int h, float pixel_size) { return ((h - 1.0f) / 2.0f - i) * pixel_size; }
float degree_to_rad(float theta_degree) {	return theta_degree * M_PI / 180.; }
float rad_to_degree(float theta) { return theta * 180. / M_PI; }

Eigen::Vector3f calculate_unit_vector(Eigen::Vector3f past, Eigen::Vector3f curr)
{
	Eigen::Vector3f d = curr - past;
	float d_norm = d.norm();
	d /= d_norm;
	return d;
}
