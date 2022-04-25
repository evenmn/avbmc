#include "constraint.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   Constraint base class
------------------------------------------------------------------------------- */

Constraint::Constraint(Box* box_in)
{
    box = box_in;
}
