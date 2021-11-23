#include "particle.h"

Particle::Particle(const std::string label_in, const std::valarray<double> r_in)
{
    /* Particle constructor
     */
    label = label_in;
    r = r_in;
}

