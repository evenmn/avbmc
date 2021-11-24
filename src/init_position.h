#pragma once
#include <iostream>
#include <vector>
#include <valarray>
#include <string>

using namespace std;

std::vector<std::valarray<double> > fcc(const int ncells, const double lenbulk, const int ndim=3);
std::vector<class Particle *> from_xyz(const std::string filename);
