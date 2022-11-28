#include <iostream>
#include <algorithm>

#include "rng.h"

RandomNumberGenerator::RandomNumberGenerator() {}

/* ----------------------------------------------------------------------------
   Inspired by numpy.random.choice, but simplified as we always want to draw
   only one sample. The samples are also always range(len(probabilities), so we
   don't need to take them as an argument.
---------------------------------------------------------------------------- */
/*
int choice(const std::vector<double> &probabilities)
{
    double r = next_double();
    double p = 0.;
    for (unsigned int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}
*/

/* ----------------------------------------------------------------------------
   Shuffle vector randomly
---------------------------------------------------------------------------- */

std::vector<unsigned int> RandomNumberGenerator::shuffle(std::vector<unsigned int> &vector)
{
    std::random_shuffle(vector.begin(), vector.end());
    return vector;
}
