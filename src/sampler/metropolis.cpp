#include "metropolis.h"
#include "../box.h"


Metropolis::Metropolis(class Box* box_in)
    : Sampler(box_in) {}

double Metropolis::get_acceptance_prob(class Moves* move, double temp, double chempot)
{
    /* Get acceptance probability of Metropolis sampling
     */
    return move->accept() * exp((du + chempot)/temp);
}
