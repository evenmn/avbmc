#include <iostream>

#include "boundary.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   Boundary base class constructor, overwriting box pointer
------------------------------------------------------------------------------- */

Boundary::Boundary(Box* box_in)
{
    box = box_in;
}


/* ----------------------------------------------------------------------------
   If the correct_position is not overwritten by subclass, the default
   behaviour is to return true.
------------------------------------------------------------------------------- */

bool Boundary::correct_position()
{
    return true;
}


/* ----------------------------------------------------------------------------
   If the correct_velocity is not overwritten by subclass, the default
   behaviour is to return true.
------------------------------------------------------------------------------- */

bool Boundary::correct_velocity()
{
    return true;
}

/* ----------------------------------------------------------------------------
   If the correct_distance is not overwritten by subclass, the default
   behaviour is to return true.
------------------------------------------------------------------------------- */

bool Boundary::correct_distance()
{
    return true;
}
