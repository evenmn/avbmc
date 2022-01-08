#pragma once
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

using namespace std;


class Thermo
{
public:
    Thermo(class Box* box_in, const int freq_in, const std::string filename, const std::vector<std::string> outputs_in);
    void print_header();
    void print_line(int);
    ~Thermo();

private:
    class Box* box = nullptr;

    std::vector<std::function<double(class Box*)> > output_functions;
    std::vector<std::string> outputs;

    int freq;

    std::ofstream f;
};
