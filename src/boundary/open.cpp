#include <iostream>

#include "../box.h"
#include "open.h"


/* ----------------------------------------------------------------------------
   Open boundary constructor
------------------------------------------------------------------------------- */

Open::Open(Box* box_in)
    : Boundary(box_in)
{
    label = "Open";
}
