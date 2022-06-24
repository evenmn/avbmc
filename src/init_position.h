/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <vector>
#include <valarray>
#include <string>


std::vector<std::valarray<double> > fcc(unsigned int, double, unsigned int = 3);
std::vector<class Particle *> from_xyz(const std::string &);
