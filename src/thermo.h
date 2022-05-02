#pragma once
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <vector>


class Thermo
{
public:
    Thermo(class Box *, int, const std::string &, const std::vector<std::string> &);
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
