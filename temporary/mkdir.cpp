#include "fileio.h"
// #include <stdlib.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <vector>


int main(){
	float i = 1.5;

   // std::filesystem::create_directories("sandbox/1/2/a");
   // std::filesystem::create_directory("sandbox/1/2/b");
	// std::string hoge = "hoge";
	// mkdir(hoge.c_str(), 0777);
	// mkdir("hoge/test", 0777);

	std::string path;
	std::ostringstream ostr;
	ostr << "result_" << "Shepp" << "/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
	std::cout << path.c_str() << std::endl;
	ostr << "from_" << "origin" << "/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
	std::cout << path.c_str() << std::endl;
	ostr << i << "mm_collimator/";
	path = ostr.str();
	mkdir(path.c_str(), 0777);
	std::cout << path.c_str() << std::endl;

}