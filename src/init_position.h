#pragma once
#include <iostream>
#include <vector>
#include <valarray>
#include <string>


std::vector<std::valarray<double> > fcc(int, double, int = 3);
std::vector<class Particle *> from_xyz(const std::string &);
