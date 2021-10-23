#include "euler.h"
#include "../box.h"


Euler::Euler(class Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}

void Euler::next_step()
{
    /* Move particle one step
     */

    box->positions += box->velocities * dt;
    box->velocities += box->accelerations * dt;
    box->forcefield->eval_acc(box->positions, box->accelerations);
}
