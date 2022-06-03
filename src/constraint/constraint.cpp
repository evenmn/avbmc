/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

#include "constraint.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   Constraint base class
------------------------------------------------------------------------------- */

Constraint::Constraint(Box* box_in)
{
    box = box_in;
    cutoff_id = type1 = type2 = 0;
}


/* ----------------------------------------------------------------------------
   Copy constructor
------------------------------------------------------------------------------- */

Constraint::Constraint(const Constraint &other) : neigh_list(other.neigh_list)
{
    cutoff_id = other.cutoff_id;
    type1 = other.type1;
    type2 = other.type2;
    box = other.box;
}
