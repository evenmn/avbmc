#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>


using namespace std;

std::vector<std::string> split(const std::string s);
std::vector<class Particle *> read_xyz(const std::string filename);
void write_xyz(std::ofstream&, double **, const int, const int, const std::vector<std::string>, const std::string);

