#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"

#define W 2048
#define H 1024
#define COUNT 8

struct Condition {
	int detector_num;
	int photon_num;
	int detector_size_w;
	int detector_size_h;
	int surface_source_w;
	int surface_source_h;
	float rotation_radius;
	float distance_collimator_to_detector;
	float pixel_size_detector;
	float collimator_h;
	float collimator_w;
	int collimator_shape_int; // 1: circle , 2: square
};

int isFOV()
{
	int i = 0;

	return i ? 1 : 0;
}


int main()
{
	Condition cond;
	
	cond.collimator_shape_int = 1; // 1: circle , 2: square

	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.photon_num = 100000000;
	// cond.photon_num = 1000000000;
	cond.surface_source_w = W;
	cond.surface_source_h = H;
	cond.distance_collimator_to_detector = 7.5;
	cond.pixel_size_detector = 0.08;
	cond.collimator_h = 1.;
	cond.collimator_w = 0.5;
	cond.rotation_radius = 25.;


	
	float x1 = cond.collimator_w / 2.;
	float y1 = -cond.collimator_h / 2.;
	float x2 = cond.collimator_w / 2. + (cond.collimator_h / 2.) * tan(M_PI / 6.);
	float y2 = -cond.collimator_h;
	float t = - x1 / (x2 - x1);
	float y = y1 + t * (y2 - y1);

	Eigen::Vector3f cross_point_edge;
	cross_point_edge << 0, y, 0;

	
	std::cout << "点の位置は" << std::endl;
	std::cout << "y = " << y << std::endl;

	std::cout << "t * (y2 - y1) = " << t * (y2 - y1) << std::endl;

	/*---- 有効視野の計算 ----*/
	// Eigen::Vector3f d_vec;
	// d_vec << 0, -1., 0;
	// Eigen::Vector3f max_vec;
	// max_vec << (cond.collimator_w / 2. + (cond.collimator_h / 2.) * tan(M_PI / 6.)), -cross_point_edge(1) + (-(cond.collimator_h / 2. + cond.rotation_radius)) , 0.;
	// float max_cos_theta = d_vec.dot(max_vec) / max_vec.norm();
	// float max_theta = acos(max_cos_theta);
	/*---------------------*/

	// std::cout << "d_vec = " << d_vec << std::endl;
	// std::cout << "max_vec = " << max_vec << std::endl;
	// std::cout << "max_vec.norm() = " << max_vec.norm() << std::endl;
	// std::cout << "max_cos_theta = " << max_cos_theta << std::endl;


	// std::cout << "max_theta = " << max_theta << std::endl;
	// std::cout << "max_theta_degree = " << max_theta * 180 / M_PI << std::endl;
	// std::cout << "next\n\n" <<  std::endl;


	
	Eigen::Vector3f on_detector;
	on_detector << 8, -( cond.rotation_radius + cond.distance_collimator_to_detector), 0.;

	
	float max_theta = M_PI / 6.;
	Eigen::Vector3f p_vec;
	p_vec << on_detector(0) - cross_point_edge(0), on_detector(1) - cross_point_edge(1), on_detector(2) - cross_point_edge(2); 
	Eigen::Vector3f d_vec;
	d_vec << 0, -1., 0;

	float cos_theta = d_vec.dot(p_vec) / p_vec.norm();
	float theta = acos(cos_theta);

	std::cout << "theta = " << theta << std::endl;
	int isfov = theta < max_theta ? 1 : 0;
	std::cout << "judge = " << isfov << std::endl;
	



}