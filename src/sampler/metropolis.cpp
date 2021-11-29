#include "metropolis.h"
#include "../box.h"
#include "../moves/moves.h"


Metropolis::Metropolis(Box* box_in)
    : Sampler(box_in) {}

double Metropolis::w(const int npar)
{
    /* Get weight function evaluation
     */
    return 1.;
}
