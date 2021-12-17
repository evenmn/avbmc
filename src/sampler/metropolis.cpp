#include "metropolis.h"
#include "../box.h"
#include "../moves/moves.h"


/* ---------------------------------------
   Metropolis constructor
------------------------------------------ */

Metropolis::Metropolis(Box* box_in)
    : Sampler(box_in) {}


/* ---------------------------------------
   Weight function (umbrella) to be added
   to energy for biased moves. For 
   Metropolis, this should always be zero
------------------------------------------ */

double Metropolis::w(const int npar)
{
    return 0.;
}
