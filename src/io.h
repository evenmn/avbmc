#pragma once
#include <fstream>
#include <string>
#include <vector>


std::vector<std::string> split(std::string);
std::vector<class Particle *> read_xyz(std::string);
void write_xyz(std::ofstream&, double **, int, int, std::vector<std::string>, std::string);

