#include "metropolis.h"
#include "../system.h"


/* ---------------------------------------
   Metropolis constructor
------------------------------------------ */

Metropolis::Metropolis(System* system_in)
    : Sampler(system_in) {}


/* ---------------------------------------
   Weight function (umbrella) to be added
   to energy for biased moves. For 
   Metropolis, this should always be zero
------------------------------------------ */

double Metropolis::w(const int npar)
{
    return 0.;
}
