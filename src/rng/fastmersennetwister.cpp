/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Faster implementation of Mersenne Twister, see
  https://github.com/cslarsen/mersenne-twister
---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <ctime>

#include "fastmersennetwister.h"

namespace mt {
    #include "cslarsen-mersenne-twister.h"
}


/* ----------------------------------------------------------------------------
   FastMersenneTwister constructor
---------------------------------------------------------------------------- */

FastMersenneTwister::FastMersenneTwister(int seed_) : RandomNumberGenerator()
{
    label = "Fast-Mersenne-Twister";

    if (seed_ > 0) {
        mt::seed(seed_);
    }
    else {
        mt::seed(std::time(NULL));   
    }
}


/* ----------------------------------------------------------------------------
   Set seed
---------------------------------------------------------------------------- */

void FastMersenneTwister::set_seed(unsigned int seed)
{
    mt::seed(seed);
}


/* ----------------------------------------------------------------------------
   Returns a double drawn from a Gaussian distribution with mean 'mean' and
   variance 'variance'
---------------------------------------------------------------------------- */

double FastMersenneTwister::next_gaussian(double mean, double variance)
{
    std::cout << "No gaussian numbers" << std::endl;
    exit(0);
}


/* ----------------------------------------------------------------------------
   Returns an integer between 0 and (including) 'upper_limit'
---------------------------------------------------------------------------- */

int FastMersenneTwister::next_int(int upper_limit)
{
    return mt::rand_u32() % upper_limit;
}

uint32_t FastMersenneTwister::next_int(const uint32_t &upper_limit)
{
    return (mt::rand_u32() % upper_limit);
}


/* ----------------------------------------------------------------------------
   Returns a double between 0 and 1 drawn from a uniform distribution
---------------------------------------------------------------------------- */

double FastMersenneTwister::next_double()
{
    return static_cast <double> (mt::rand_u32()) * RAND_MAX_inv;
    //return (mt::rand_u32() * RAND_MAX_inv);
}


/* ----------------------------------------------------------------------------
   Inspired by numpy.random.choice, but simplified as we always want to draw
   only one sample. The samples are also always range(len(probabilities), so we
   don't need to take them as an argument.
---------------------------------------------------------------------------- */

int FastMersenneTwister::choice(const std::vector<double> &probabilities)
{
    double r = next_double();
    double p = 0.;
    for (unsigned int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}
