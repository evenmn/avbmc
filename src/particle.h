/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>
#include <valarray>


class Particle
{
public:
    Particle(const std::string &, const std::valarray<double> &);
    Particle(const Particle &);

    std::valarray<double> r;
    std::valarray<double> v;
    std::valarray<double> f;

    int type;
    double mass, poteng;
    std::string label;
};
