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
