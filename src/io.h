/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <fstream>
#include <string>
#include <vector>


std::vector<std::string> split(const std::string &);
std::vector<class Particle> read_xyz(const std::string &);
void write_xyz(std::ofstream&, double **, unsigned int, unsigned int, const std::vector<std::string> &, const std::string &);
void write_vector(std::vector<int>, const std::string &, const std::string &);
void write_vector(std::vector<double>, const std::string &, const std::string &);
void write_array(int *, int length, const std::string &, const std::string &);
void write_array(unsigned int *, unsigned int length, const std::string &, const std::string &);
