/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Mersenne-Twister pseudo random number generator (prng). The prng has a 
  period of 2^19937, and is popularly abbreviated mt19937. The algorithm was
  first proposed and described by ... This code is just a wrapper around the
  standard C++ implementation of Mersenne-Twister.
---------------------------------------------------------------------------- */

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>

#include "mersennetwister.h"


/* ----------------------------------------------------------------------------
   Mersenne Twister constructor
---------------------------------------------------------------------------- */

MersenneTwister::MersenneTwister(int seed_)
    : RandomNumberGenerator(), generator(seed())
{
    label = "Mersenne-Twister";
    if (seed_ > 0) {
        generator.seed (seed_);
    }

    //std::uniform_real_distribution randu_(0, 1);
    //randu = randu_;
}


/* ----------------------------------------------------------------------------
   Set seed
---------------------------------------------------------------------------- */

void MersenneTwister::set_seed(unsigned int seed)
{
    generator.seed (seed);
}


/* ----------------------------------------------------------------------------
   Returns a double drawn from a Gaussian distribution with mean 'mean' and
   variance 'variance'
---------------------------------------------------------------------------- */

double MersenneTwister::next_gaussian(double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return dis(generator);
}


/* ----------------------------------------------------------------------------
   Returns an interger between 0 and (including) 'upper_limit'
---------------------------------------------------------------------------- */

int MersenneTwister::next_int(int upper_limit)
{
    std::uniform_int_distribution<int> dis(0, upper_limit - 1);
    return dis(generator);
}


/* ----------------------------------------------------------------------------
   Returns a double between 0 and 1 drawn from a uniform distribution
---------------------------------------------------------------------------- */

double MersenneTwister::next_double()
{
    //std::uniform_real_distribution<double> dis(0, 1);
    //return dis(generator);
    return randu(generator);
}


/* ----------------------------------------------------------------------------
   Inspired by numpy.random.choice, but simplified as we always want to draw
   only one sample. The samples are also always range(len(probabilities), so we
   don't need to take them as an argument.
---------------------------------------------------------------------------- */

int MersenneTwister::choice(const std::vector<double> &probabilities)
{
    double r = next_double();
    double p = 0.;
    for (unsigned int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}


/* ----------------------------------------------------------------------------
   Shuffle vector randomly
---------------------------------------------------------------------------- */

std::vector<unsigned int> MersenneTwister::shuffle(std::vector<unsigned int> &vector)
{
    std::shuffle(vector.begin(), vector.end(), generator);
    return vector;
}
