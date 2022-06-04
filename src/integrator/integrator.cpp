/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  The integrators are used to go from one molecular dynamics state to another.
---------------------------------------------------------------------------- */

#include <iostream>

#include "integrator.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   Integrator base class constructor
---------------------------------------------------------------------------- */

Integrator::Integrator(Box* box_in, double dt_in)
{
    box = box_in;
    dt = dt_in;
}
