#pragma once
#include <string>
#include <vector>


void parser(int, char**);
void set(class System &, class Box &, std::vector<std::string>, int);
void add(class System &, class Box &, std::vector<std::string>, int);
void take(class System &, class Box &, std::vector<std::string>, int);
void run(class System &, class Box &, std::vector<std::string>, int);
void thermo(class System &, class Box &, std::vector<std::string>, int);
void dump(class System &, class Box &, std::vector<std::string>, int);
void write(class System &, class Box &, std::vector<std::string>, int);
void rm(class System &, class Box &, std::vector<std::string>, int);
