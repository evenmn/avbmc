#include <iostream>

#include "integrator.h"
#include "../box.h"


/* ----------
   Integrator base class
------- */

Integrator::Integrator(Box* box_in, double dt_in)
{
    box = box_in;
    dt = dt_in;
}
