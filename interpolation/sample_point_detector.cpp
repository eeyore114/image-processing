
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include "fileio.h"

typedef struct {
	int img_w;
	int img_h;
	int img_d;
	int detector_num;
	int detector_size_w;
	int detector_size_h;
	int update_count;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
	float img_pixel_size;
	float d_width;
	float d_height;
	float time;
} Condition;

void sample_point_interpolation(std::vector<float> &detector, std::vector<float> &spi_detector, Condition &cond);

int main()
{
  Condition cond;
	cond.img_w = 128;
	cond.img_h = 128;
	cond.img_d = 128;
	cond.detector_num = 180;
	cond.detector_size_w = 512;
	cond.detector_size_h = 256;
	cond.rotation_radius = 12;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.2;
	cond.update_count = 10;
	cond.img_pixel_size = 0.1;
	cond.d_width = 0.08;
	cond.d_height = 0.08;

  std::vector<float> detector(cond.detector_num * cond.detector_size_w * cond.detector_size_h);
  std::vector<float> spi_detector(cond.detector_num * cond.detector_size_w * 2 * cond.detector_size_h * 2);
  readRawFile("efficiency_correction_float_512-256-180.raw", detector);
  sample_point_interpolation(detector, spi_detector, cond);
  writeRawFile("eefficiency_correction_float_1024-512-180.raw", spi_detector);
  std::cout << "cond.detector_size_h = " << cond.detector_size_h << std::endl;

}


void sample_point_interpolation(std::vector<float> &detector, std::vector<float> &spi_detector, Condition &cond)
{
  for(int theta_degree = 0; theta_degree < cond.detector_num; theta_degree++)
  {
    for(int i = 0; i < cond.detector_size_h; i++)
  	{
  		for(int j = 0; j < cond.detector_size_w; j++)
  		{
  			int spi_detector_w = cond.detector_size_w * 2;
  			int spi_detector_h = cond.detector_size_h * 2;
  			int spi_i = i * 2;
  			int spi_j = j * 2;
  			spi_detector[theta_degree * spi_detector_w * spi_detector_h + spi_i * spi_detector_w + spi_j] 					  = detector[theta_degree * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w + j];
  			spi_detector[theta_degree * spi_detector_w * spi_detector_h + spi_i * spi_detector_w + (spi_j + 1)]			  = detector[theta_degree * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w + j];
  			spi_detector[theta_degree * spi_detector_w * spi_detector_h + (spi_i + 1) * spi_detector_w + spi_j]			  = detector[theta_degree * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w + j];
  			spi_detector[theta_degree * spi_detector_w * spi_detector_h + (spi_i + 1) * spi_detector_w + (spi_j + 1)] = detector[theta_degree * cond.detector_size_w * cond.detector_size_h + i * cond.detector_size_w + j];
  	  }
    }
  }
  cond.detector_size_h *= 2;
}
