#include <iostream>
#include <string>
#include <valarray>

#include "particle.h"

/* ----------------------------------------------------------------------------
   Particle constructor. Takes particle label 'label_in' and particle initial
   position 'r_in' as arguments.
-------------------------------------------------------------------------------
*/

Particle::Particle(const std::string label_in, const std::valarray<double> r_in)
{
    label = label_in;
    r = r_in;
}

