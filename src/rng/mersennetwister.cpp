#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>

#include "mersennetwister.h"


/* ---------------------------------------------------------
   Mersenne Twister constructor
------------------------------------------------------------ */

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator(), generator(seed())
{
    label = "Mersenne-Twister";
}


/* ---------------------------------------------------------
   Returns a double drawn from a Gaussian distribution with 
   mean 'mean' and variance 'variance'
------------------------------------------------------------ */

double MersenneTwister::next_gaussian(double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return dis(generator);
}


/* ---------------------------------------------------------
   Returns an interger between 0 and (including) 'upper_limit'
------------------------------------------------------------ */

int MersenneTwister::next_int(int upper_limit)
{
    std::uniform_int_distribution<int> dis(0, upper_limit - 1);
    return dis(generator);
}


/* ---------------------------------------------------------
   Returns a double between 0 and 1 drawn from a uniform 
   distribution
------------------------------------------------------------ */

double MersenneTwister::next_double()
{
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}


/* ---------------------------------------------------------
   Inspired by numpy.random.choice, but simplified as we always
   want to draw only one sample. The samples are also always 
   range(len(probabilities), so we don't need to take them as
   an argument
------------------------------------------------------------ */

int MersenneTwister::choice(std::vector<double> probabilities)
{
    unsigned int i;
    double r, p;

    r = next_double();
    p = 0.;
    for (i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}
