/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Simple C++ random numbers. This should only be used to benchmark the
  performance of other pseudo random number generators, as it might be
  biased
---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>

#include "simple.h"


/* ----------------------------------------------------------------------------
   Simple constructor
---------------------------------------------------------------------------- */

Simple::Simple(int seed_) : RandomNumberGenerator()
{
    label = "Simple";

    std::srand (static_cast <unsigned> (time(NULL)));
    if (seed_ > 0) {
        std::srand (seed_);
    }
}


/* ----------------------------------------------------------------------------
   Set seed
---------------------------------------------------------------------------- */

void Simple::set_seed(unsigned int seed)
{
    std::srand  (seed);
}


/* ----------------------------------------------------------------------------
   Returns a double drawn from a Gaussian distribution with mean 'mean' and
   variance 'variance'
---------------------------------------------------------------------------- */

double Simple::next_gaussian(double mean, double variance)
{
    std::cout << "No gaussian numbers" << std::endl;
    exit(0);
}


/* ----------------------------------------------------------------------------
   Returns an interger between 0 and (including) 'upper_limit'
---------------------------------------------------------------------------- */

int Simple::next_int(int upper_limit)
{
    assert (upper_limit < RAND_MAX);
    return rand() % upper_limit;
}


/* ----------------------------------------------------------------------------
   Returns a double between 0 and 1 drawn from a uniform distribution
---------------------------------------------------------------------------- */

double Simple::next_double()
{
    return static_cast <double> (rand()) * RAND_MAX_inv;
}


/* ----------------------------------------------------------------------------
   Inspired by numpy.random.choice, but simplified as we always want to draw
   only one sample. The samples are also always range(len(probabilities), so we
   don't need to take them as an argument.
---------------------------------------------------------------------------- */

int Simple::choice(const std::vector<double> &probabilities)
{
    double r = next_double();
    double p = 0.;
    for (unsigned int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}
