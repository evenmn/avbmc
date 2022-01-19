/* ----------------------------------------------------------------------------
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Metropolis sampling, as originally proposed by Metropolis, 
   Rosenbluth, Rosenbluth, Teller and Teller (1954). The acceptance
   of the move is decided by the move's 'accept'-method, as the
   acceptance criterion is move-specific
------------------------------------------------------------------------------- */

#include "metropolis.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   Metropolis constructor
------------------------------------------------------------------------------- */

Metropolis::Metropolis(System* system_in)
    : Sampler(system_in) 
{
    label = "Metropolis";
}


/* ----------------------------------------------------------------------------
   Weight function (umbrella) to be added
   to energy for biased moves. For 
   Metropolis, this should always be zero
------------------------------------------------------------------------------- */

double Metropolis::w(const int /*npar*/)
{
    return 0.;
}
