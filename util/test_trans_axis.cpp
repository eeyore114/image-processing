// なんで座標変換の時に-1するの?


#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

int main()
{
	int i, j, w, h;
	float x, y;
	float pixel_size = 1.;
	w = h = 129;


	float test = 1.0;

	x = test;
	y = 0;
	j = floor( x / pixel_size + (w - 1.0f) / 2.0f );
	i = floor( (h - 1.0f) / 2.0f - y / pixel_size );

	x = (j - (w) / 2.0f) * pixel_size;
	y = ((h) / 2.0f - i) * pixel_size;
	// x = (j - (w - 1.0f) / 2.0f) * pixel_size;
	// y = ((h - 1.0f) / 2.0f - i) * pixel_size;

	std::cout << "------- -1した -------" << '\n';
	std::cout << "j = " << j << '\n';
	std::cout << "i = " << i << '\n';
	std::cout << "x = " << x << '\n';
	std::cout << "y = " << y << '\n';

	x = test;
	y = 0;

	j = floor( x / pixel_size + (w) / 2.0f );
	i = floor( (h) / 2.0f - y / pixel_size );
	x = (j - (w) / 2.0f) * pixel_size;
	y = ((h) / 2.0f - i) * pixel_size;
	// x = (j - (w - 1.0f) / 2.0f) * pixel_size;
	// y = ((h - 1.0f) / 2.0f - i) * pixel_size;

	std::cout << "------- してない -------" << '\n';
	std::cout << "j = " << j << '\n';
	std::cout << "i = " << i << '\n';
	std::cout << "x = " << x << '\n';
	std::cout << "y = " << y << '\n';
}
