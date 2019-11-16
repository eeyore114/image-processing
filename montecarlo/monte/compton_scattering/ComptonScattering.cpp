#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#define N 4


void ComptonScattering(float &scattering_energy, float &scattering_angle);


template <class T>
T mkrandom(T min, T max);

int main()
{
	int count = 0;
	for(int i = 0; i < N; i++)
	{
		float scattering_energy;
		float scattering_angle;

		ComptonScattering(scattering_energy, scattering_angle);
		if(scattering_energy != 0. && scattering_angle != 0.)
		{
			cout << "energy = " << scattering_energy << " angle = " << scattering_angle << endl;
			count++;
		}
	}
	cout << "count = " << count << endl;
}

void ComptonScattering(float &scattering_energy, float &scattering_angle)
{
	float h = 6.62 * pow(10., -34.);
	float m0 = 9.11 * pow(10., -28.);
	float c = 3.0 * pow(10., 8.);
	float E = 140.0;

	float lambda = m0 * pow(c, 2.) / E;
	float relative_frequency = (lambda + 2.) / (9 * lambda + 2);
	// cout << m0 + pow(c , 2.) / 3.639 << "\n";
	
	cout << "lambda = " << lambda << "  relative_frequency = " << relative_frequency << endl;
	
	float r1 = mkrandom(0., 1.);
	float r2 = mkrandom(0., 1.);
	float r3 = mkrandom(0., 1.);
	
	float rho, cmp;
	if(relative_frequency < r1)
	{
		rho = 1. + (2. / lambda) * r2;
		cmp = 4 * (1. / rho - 1 / pow(rho, 2.));
	}
	else
	{
		rho = (2. + lambda) / (lambda + 2. * (1. - r2));
		cmp = (pow(lambda - rho * lambda + 1., 2.) + 1 / rho) / 2.;
	}

	if(r3 <= cmp)
	{
		float lambda_dash = rho * lambda;
		float tmp1 = 1. - m0 * c * (lambda_dash - lambda) / h;
		float tmp2 = 1. - (lambda_dash - lambda);
		scattering_energy = m0 * pow(c, 2.) / lambda_dash;
		float theta = acos(tmp2);
		scattering_angle = theta * 180.0f / M_PI;

		// cout << "lambda_dash = " << lambda_dash << "  energy = " << scattering_energy << " angle = " << scattering_angle << endl;
	}
	else
	{
		// 棄却処理
		scattering_energy = 0.;
		scattering_angle = 0.;
	}
}




template <class T>
T mkrandom(T min, T max)
{
	random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値

    if(typeid(int) == typeid(T))
    {
        uniform_int_distribution<T> random(min, max); 
        T r = random(mt);
        return r;
    }
    else
    {
        uniform_real_distribution<T> random(min, max); 
        T r = random(mt);
        return r;
    }

}
