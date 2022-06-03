/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <random>
#include <vector>
#include <string>


class RandomNumberGenerator
{
public:
    RandomNumberGenerator() {};
    virtual void set_seed(unsigned int) = 0;
    virtual int next_int(int) = 0;
    virtual double next_double() = 0;
    virtual double next_gaussian(double = 0.0, double = 1.0) = 0;
    virtual int choice(std::vector<double>) = 0;
    virtual std::vector<unsigned int> shuffle(std::vector<unsigned int> &) = 0;
    virtual ~RandomNumberGenerator() = default;

    std::string label;

protected:
    std::random_device seed;
};
