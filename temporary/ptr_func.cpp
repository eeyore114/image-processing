#include <iostream>
using namespace std;
#include "fileio.h"

int add(int a, int b){
    return a+b;
}

int sub(int a, int b){
    return a-b;
}


int main(){

   int (*const FP_AB[])(int a, int b) = {
        add, sub
    };
    cout << FP_AB[0](1, 1) << endl;
    cout << FP_AB[1](1, 1) << endl;

}