#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;


class Thermo
{
public:
    Thermo(class Box* box_in, const int freq_in, const string filename, const vector<string> outputs_in);
    void print_line();
    void print_header();

private:
    class Box* box = nullptr;

    vector<function<double(class Box*)>> output_functions;
    vector<string> outputs;

    int freq;
};
