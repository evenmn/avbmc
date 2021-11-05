#include "eulercromer.h"
#include "../box.h"


EulerCromer::EulerCromer(class Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}

void EulerCromer::next_step()
{
    /* Move particle one step
     */

    box->velocities += box->accelerations * dt;
    box->positions += box->velocities * dt;
    box->forcefield->eval_acc(box->positions, box->accelerations, box->potengs, true);
}
