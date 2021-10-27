#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;


class Dump
{
public:
    Dump(class Box* box_in, const int freq_in, const string filename, const vector<string> outputs_in);
    void print_frame();
    ~Dump();

private:
    class Box* box = nullptr;

    vector<function<mat(class Box*)>> output_functions;
    vector<string> outputs;

    int freq;
    ofstream f;
    string info_line;
};
