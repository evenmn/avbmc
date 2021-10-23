#include "integrator.h"
//#include "../box.h"


Integrator::Integrator(class Box* box_in, double dt_in)
{
    box = box_in;
    dt = dt_in;
}

//Integrator::~Integrator() {}
