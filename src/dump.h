#pragma once
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <vector>


class Dump
{
public:
    Dump(class Box*, int, std::string, std::vector<std::string>);
    void print_frame(int);
    ~Dump();

private:
    class Box* box = nullptr;

    std::vector<std::function<double **(class Box *)> > output_functions;
    std::vector<std::string> outputs;
    std::vector<int> nvars;

    int freq, nvar;
    std::ofstream f;
    std::string info_line;
};
