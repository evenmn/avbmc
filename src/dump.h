/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <vector>


class Dump
{
public:
    Dump(class Box *, unsigned int, const std::string &, std::vector<std::string>);
    void print_frame(unsigned int);
    ~Dump();

private:
    class Box* box = nullptr;

    std::vector<std::function<double **(class Box *)> > output_functions;
    std::vector<std::string> outputs;
    std::vector<int> nvars;

    unsigned int freq, nvar;
    std::ofstream f;
    std::string info_line;
};
