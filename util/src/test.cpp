/*

--------------
	util概要
--------------

・汎用的に使える関数をまとめておくもの
・必ずテストコードを書く
・cudaで使う場合は別途でutil_deviceを作成



*/

#include <stdio.h>
#include <unistd.h>
#include <iostream>
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


bool test_is_out_of_image();
bool test_calculate_unit_vector();
void output_test_result(int success, int failed);



int main()
{
	int success = 0;
	int failed = 0;

	if(test_is_out_of_image()) { success++; } else { failed++; }
	output_test_result(success, failed);
}


bool test_is_out_of_image()
{
	std::cout << "test is_out_of_image function" << '\n';
	int i, j, h, w;

	/*---- 1 ----*/
	i = 127; j = 127; h = 128; w = 128;
	if(is_out_of_image(i, h, j, w) == true)
	{
		std::cout << "\n-> test_is_out_of_image is failed. 1 \n" << '\n';
		return false;
	}

	/*---- 2 ----*/
	i = -1; j = 127; h = 128; w = 128;
	if(is_out_of_image(i, h, j, w) == false)
	{
		std::cout << "\n-> test_is_out_of_image is failed. 2 \n" << '\n';
		return false;
	}
	return true;
}

//
// bool test_calculate_unit_vector()
// {
// 	int a = 1;
// }



//
// bool test_transform_img_coordinate_same_axis()
// {
//
// }
//
// bool test_transform_img_coordinate_opposite_axis()
// {
//
// }
//
//
// bool transform_origin_coordinate_same_axis()
// {
//
// }
//
// bool test_transform_origin_coordinate_opposite_axis()
// {
//
// }




void output_test_result(int success, int failed)
{
	printf("\x1b[32m");     /* 前景色を緑に */
  printf("success : %d", success);
	printf("\x1b[39m");     /* 前景色をデフォルトに戻す */
	if(failed > 0)
	{
		printf(", ");
		printf("\x1b[31m");     /* 前景色を赤に */
    printf("failed : %d\n", failed);
	}
	else { printf("\n"); }
}
