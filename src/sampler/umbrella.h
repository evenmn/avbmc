/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <functional>
#include <valarray>

#include "sampler.h"


class Umbrella : public Sampler
{
public:
    Umbrella(class System*, std::function<double (int)>, int = 100);
    Umbrella(class System*, std::valarray<double>);
    double w(int) override;

private:
    std::valarray<double> tabulate(std::function<double (int)>, int);
    std::valarray<double> tabulated;
    int maxpar;
    std::function <double (int)> f;
};
