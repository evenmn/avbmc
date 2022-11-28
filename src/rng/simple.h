/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <vector>
#include <cstdlib>

#include "rng.h"


class Simple : public RandomNumberGenerator
{
public:
    Simple(int = -1);
    void set_seed(unsigned int) override;
    int next_int(int) override;
    double next_double() override;
    double next_gaussian(double, double) override;
    int choice(const std::vector<double> &) override;
    //std::vector<unsigned int> shuffle(std::vector<unsigned int> &) override;
    ~Simple() = default;

private:
    static constexpr double RAND_MAX_inv = 1. / static_cast <double> (RAND_MAX);
};
