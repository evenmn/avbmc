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
