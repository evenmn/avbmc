/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Metropolis sampling, as originally proposed by Metropolis, Rosenbluth,
   Rosenbluth, Teller and Teller (1954). The acceptance of the move is decided
   by the move's 'accept'-method, as the acceptance criterion is move-specific
---------------------------------------------------------------------------- */

#include "metropolis.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   Metropolis constructor
---------------------------------------------------------------------------- */

Metropolis::Metropolis(System* system_in)
    : Sampler(system_in) 
{
    label = "Metropolis";
}


/* ----------------------------------------------------------------------------
   Weight function (umbrella) to be added to energy for biased moves. For 
   Metropolis, this should always be zero
---------------------------------------------------------------------------- */

double Metropolis::w(const int /*npar*/)
{
    return 0.;
}
