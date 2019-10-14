#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>


template <typename T> void readRawFile(std::string fname, std::vector<T>& image);
template <typename T> void writeRawFile(std::string fname, std::vector<T>& image);

template <typename T> void writeRawFile(std::string fname, std::vector<T>& image)
{
    std::ofstream ofs(fname.c_str()); 
    if (!ofs) std::cerr << "Failed to open " << fname << std::endl;
    else { ofs.write(reinterpret_cast<char*>(&image[0]), sizeof(T) * image.size()); }
}


template <typename T> void readRawFile(std::string fname, std::vector<T>& image)
{
    std::ifstream ifs(fname.c_str());
    if (!ifs) std::cerr << "Failed to open " << fname << std::endl;
    else { ifs.read(reinterpret_cast<char*>(&image[0]), sizeof(T) * image.size()); }
}
