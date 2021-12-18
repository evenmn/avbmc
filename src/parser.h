#pragma once
#include <string>
#include <vector>


void parser(int, char**);
void set(class Box &, std::vector<std::string>, int);
void add(class Box &, std::vector<std::string>, int);
void take(class Box &, std::vector<std::string>, int);
void run(class Box &, std::vector<std::string>, int);
void thermo(class Box &, std::vector<std::string>, int);
void dump(class Box &, std::vector<std::string>, int);
void write(class Box &, std::vector<std::string>, int);
void rm(class Box &, std::vector<std::string>, int);
