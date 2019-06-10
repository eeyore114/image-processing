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
#include "xcom.h"


void ComptonScattering(double* energy_counter);


int main()
{
	double* energy_counter = (double*)calloc(180, sizeof(double));

	for(int i = 0; i < 100000000; i++)
	{
    	if(i % 1000000 == 0){ printf("photon : %09d\n\n", i); }
		ComptonScattering(energy_counter);
	}

	FILE* fp;
    if((fp = fopen("compton_test.csv","w")) != NULL)
    {
        for(int i = 0; i < 180; i++)
        {
	        fprintf(fp,"%d,%f\n", i, energy_counter[i]);
        }
        fclose(fp);
    }



}



void ComptonScattering(double* energy_counter)
{
	float h = 6.62 * pow(10., -34.);
	float m0 = 9.11 * pow(10., -31.);
	float c = 3.0 * pow(10., 8.);
	float KeV_to_J = 1.602 * pow(10., -19.) * pow(10., 3.);

	float energy = 140. * KeV_to_J;

	float lambda = m0 * pow(c, 2.) / energy;
	float relative_frequency = (lambda + 2.) / (9 * lambda + 2);
	float rho;
	int flag = 0;

	while(flag == 0)
	{
		float r1 = genrand_real1();
		float r2 = genrand_real1();
		float r3 = genrand_real1();

		if(relative_frequency < r1)
		{
			rho = 1. + (2. / lambda) * r2;
			flag = (r3 < 4 * (1. / rho - 1 / pow(rho, 2.)));
		}
		else
		{
			rho = (2. + lambda) / (lambda + 2. * (1. - r2));
			flag = (r3 < (pow(lambda - rho * lambda + 1., 2.) + 1 / rho) / 2.);
		}
	}

	float lambda_dash = rho * lambda;
	energy = m0 * pow(c, 2.) / lambda_dash;
	energy /= KeV_to_J;

	float theta_relative_cos = 1. - (lambda_dash - lambda);


	int energy_int = round(energy);
	int theta_int = round(acos(theta_relative_cos) * 180. / M_PI);
	energy_counter[theta_int] = energy_int;
}