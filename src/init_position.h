#pragma once
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat fcc(int ncells, double lenbulk, int ndim=3);
mat from_xyz(string filename);
