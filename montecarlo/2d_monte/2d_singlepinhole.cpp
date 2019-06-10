#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"

template <class T>
void readRawFile (std::string fname, const size_t num, T* image);
template <class T>
void writeRawFile (std::string fname, const size_t num, T* image);

void RaySimulation();
/*-----------Raysimulation-----------*/
void readXcom(string read_xcom_name, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o);
void WriteEnergySpectrum(float* energy_spectrum, string write_energy_spectrum_name);
void WriteImageProfile(float* detector, string write_image_profile_name, int detector_size, int detector_num);
int JudgeIsair(Eigen::Vector3f curr_, float* read_img, float scale_ratio_fantom, int mumap_size);
void WriteImage(float* detector, string medium_name, string init_position, int w_collimator_mm, string img_name, int detector_size, int detector_num);
/*-----------------------------------*/

void mlem(float* f, float* g, float* h, struct Condition cond);
/*-----------MLEM-----------*/
void make_cij(float* f, struct Condition cond);
void Projecion_SinglePinhole_3d(float* f,float* g, struct Condition cond);
void backProjection_Singlepinhole_3d(float* f, float* g, struct Condition cond);
void search_LessThan30_3d(float* f, float* g, struct Condition cond);
/*--------------------------*/




class Photon {
public:
	Eigen::Vector3f curr_;
	Eigen::Vector3f past_;
	float energy_;
	float theta_cos_, theta_sin_;
	float phi_cos_, phi_sin_;
	float optical_length_;
	float coherent_, compton_, photo_, mu_, mu_max_;
	int scatter_;
	float medium;
	int MUMAP_SIZE;

	Photon(int i, int j, int k, float scale_ratio_fantom, int mumap_size, string init_position);

	void move();
	void SetProbability(float* read_img, string medium_name, float scale_ratio_fantom, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o);
	void ComptonScattering();
};

Photon::Photon(int i, int j, int k, float scale_ratio_fantom, int mumap_size, string init_position) : scatter_(0)
{
	energy_ = 140.;
	float rnd = 1. * 2 * (genrand_real1() - 0.5);
	float theta = acos(rnd);
	theta_sin_ = sin(theta);
	theta_cos_ = cos(theta);
	float phi = genrand_real1() * 2 * M_PI;
	phi_sin_ = sin(phi);
	phi_cos_ = cos(phi);
	MUMAP_SIZE = mumap_size;

	if(init_position == "voxel")
	{
		float mu_x = - (MUMAP_SIZE - 1.0) / 2.0 + k;
		float mu_y =   (MUMAP_SIZE - 1.0) / 2.0 - j;

		mu_x /= scale_ratio_fantom;
		mu_y /= scale_ratio_fantom;

		// 初期位置にある程度ばらつきをつける
		mu_x += 0.1 * 2 * (genrand_real1() - 0.5);
		mu_y += 0.1 * 2 * (genrand_real1() - 0.5);

		past_ << mu_x, mu_y;
		curr_ << mu_x, mu_y;
	}
	else if(init_position == "origin")
	{
		past_ << 0., 0.;
		curr_ << 0., 0.;
	}
	else { cout << init_position << " can't select." << endl; exit(-1); }


}

void Photon::move()
{
	float N_per_N0 = 1. - genrand_real2();

	optical_length_ = - log(N_per_N0) / mu_max_;

	curr_(0) = past_(0) + optical_length_ * theta_sin_ * phi_cos_;
	curr_(1) = past_(1) + optical_length_ * theta_sin_ * phi_sin_;
	curr_(2) = past_(2) + optical_length_ * theta_cos_;
}

void Photon::SetProbability(float* read_img, string medium_name, float scale_ratio_fantom, float* Coherent_ca, float* Compton_ca, float* Photoelectric_ca, float* mu_ca,	float* Coherent_h2o, float* Compton_h2o, float* Photoelectric_h2o, float* mu_h2o)
{
	int index = floor(energy_);

	float ca = mu_ca[index];
	float h2o = mu_h2o[index];

	if(medium_name == "ca") { mu_max_ = ca; }
	if(medium_name == "h2o") { mu_max_ = h2o; }
	if(medium_name == "multi") { mu_max_ = ca > h2o ? ca : h2o; }
	if(medium_name == "air") { mu_max_ = 1; }

	int j = round(curr_(0) * scale_ratio_fantom + (MUMAP_SIZE - 1) / 2.0);
	int i = round((MUMAP_SIZE - 1) / 2.0 - (curr_(1) * scale_ratio_fantom));
	int k = round((MUMAP_SIZE - 1) / 2.0 - (curr_(2) * scale_ratio_fantom));

	if(k > 0 && k < MUMAP_SIZE && i > 0 && i < MUMAP_SIZE && j > 0 && j < MUMAP_SIZE)
	{
		medium = read_img[k * MUMAP_SIZE * MUMAP_SIZE + i * MUMAP_SIZE + j];

		if(medium > 0.2)
		{
			coherent_ = Coherent_ca[index];

			compton_ = Compton_ca[index];

			photo_ = Photoelectric_ca[index];

			mu_ = mu_ca[index];
		}
		else if(medium > 0.1)
		{
			coherent_ = Coherent_h2o[index];

			compton_ = Compton_h2o[index];

			photo_ = Photoelectric_h2o[index];

			mu_ = mu_h2o[index];
		}
		// 空気の場合（move()が動くように適当に設定）
		else { mu_ = 1.; }
	}
}

void Photon::ComptonScattering()
{
	float h = 6.62 * pow(10., -34.);
	float m0 = 9.11 * pow(10., -31.);
	float c = 3.0 * pow(10., 8.);
	float KeV_to_J = 1.602 * pow(10., -19.) * pow(10., 3.);

	energy_ *= KeV_to_J;

	float lambda = m0 * pow(c, 2.) / energy_;
	float relative_frequency = (lambda + 2.) / (9 * lambda + 2);
	float rho;
	int isAcceptable = 0;

	while(isAcceptable == 0)
	{
		float r1 = genrand_real1();
		float r2 = genrand_real1();
		float r3 = genrand_real1();

		if(relative_frequency > r1)
		{
			rho = 1. + (2. / lambda) * r2;
			isAcceptable = (r3 < 4 * (1. / rho - 1 / pow(rho, 2.)));
		}
		else
		{
			rho = (2. + lambda) / (lambda + 2. * (1. - r2));
			isAcceptable = (r3 < (pow(lambda - rho * lambda + 1., 2.) + 1 / rho) / 2.);
		}
	}

	float lambda_dash = rho * lambda;
	energy_ = m0 * pow(c, 2.) / lambda_dash;
	energy_ /= KeV_to_J;

	float theta_relative_cos = 1. - (lambda_dash - lambda);
	float theta_relative_sin = sqrt(1. - pow(theta_relative_cos, 2.));

	if(isnan(theta_relative_sin)) {theta_relative_sin = 0.;}

	float phi_relative = genrand_real1() * 2. * M_PI;

	float theta_cos_n = theta_cos_;
	float theta_sin_n = theta_sin_;
	float phi_cos_n = phi_cos_;
	float phi_sin_n = phi_sin_;

	theta_cos_ = -theta_sin_n * theta_relative_sin * cos(phi_relative) + theta_cos_n * theta_relative_cos;
	theta_sin_ = sqrt(1. - pow(theta_cos_, 2.));

	// 0割り対策
	if(theta_sin_ < 0.0000001)
	{
		float rnd_phi = 2 * M_PI * genrand_real2();
		phi_cos_ = cosf(rnd_phi);
		phi_sin_ = sinf(rnd_phi);
	}
	else
	{
		phi_cos_ = (theta_sin_n * phi_cos_n * theta_relative_cos + theta_cos_n * phi_cos_n * theta_relative_sin * cos(phi_relative) - phi_sin_n * theta_relative_sin * sin(phi_relative)) / theta_sin_;
		phi_sin_ = (theta_sin_n * phi_sin_n * theta_relative_cos + theta_cos_n * phi_sin_n * theta_relative_sin * cos(phi_relative) + phi_cos_n * theta_relative_sin * sin(phi_relative)) / theta_sin_;
	}
}

struct Condition {
	int detector_num;
	int detector_size;
	int img_w;
	int img_h;
	int photon_num;
	float rotation_radius;
	float distance_collimator_to_detector;
	float height_collimator;
	float width_collimator;
};

int main()
{
	/*----------------条件----------------*/
	// ----------------------
	// montecarlo simulation
	// ----------------------
	Condition cond;
	cond.img_w = 32;
	cond.img_h = 32;
	cond.detector_num = 180;
	cond.detector_size = 50;
	cond.rotation_radius = 10.;
	cond.distance_collimator_to_detector = 7.5;
	cond.height_collimator = 1.;
	cond.width_collimator = 0.5;
	cond.photon_num = 1E5;
	std::string mumap_name = "";



	// -------
	// ML-EM
	// -------
	int update_count = 20;
	std::string reconstruction_img_name = "";

	/*-----------------------------------*/



	int scatter = 6;
	float* detector = (float*)calloc(scatter * cond.detector_num * cond.detector_size, sizeof(float));
	float* proj_data = (float*)calloc(cond.detector_num * cond.detector_size, sizeof(float));
	float* reconstruction_img = (float*)calloc(cond.img_w * cond.img_h, sizeof(float));

	Raysimulation(detector, cond);
	getProjectionData(detector, proj_data, cond);
	mlem();

	writeRawFile();
}



/*-----------Raysimulation-----------*/

/*-----------------------------------*/


/*-----------MLEM-----------*/

/*--------------------------*/









template <class T>
void readRawFile (std::string fname, const size_t num, T* image)
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
void writeRawFile(std::string fname, const size_t num, T* image)
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
