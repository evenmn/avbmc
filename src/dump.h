#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


class Dump
{
public:
    Dump(class Box* box_in, const int freq_in, const std::string filename, const std::vector<std::string> outputs_in);
    void print_frame();
    ~Dump();

private:
    class Box* box = nullptr;

    std::vector<std::function<double **(class Box*)> > output_functions;
    std::vector<std::string> outputs;
    std::vector<int> nvars;

    int freq, nvar;
    std::ofstream f;
    std::string info_line;
};
