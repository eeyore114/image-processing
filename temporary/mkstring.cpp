#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main()
{
    string str;
    str.push_back('t');
    int a = 5;
    ostringstream ss;
    ss << a;
    str = str + ss.str();
    std::cout << str << std::endl;
}
