#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

mat fcc(const int ncells, const double lenbulk, const int ndim=3);
mat from_xyz(const string filename, vector<string> &chem_symbols);
