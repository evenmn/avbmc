/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <vector>
#include <random>

#include "rng.h"


class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister(int = -1);
    void set_seed(unsigned int) override;
    int next_int(int) override;
    double next_double() override;
    double next_gaussian(double, double) override;
    int choice(const std::vector<double> &) override;
    std::vector<unsigned int> shuffle(std::vector<unsigned int> &);
    ~MersenneTwister() = default;

private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> randu; //(0, 1);
};
