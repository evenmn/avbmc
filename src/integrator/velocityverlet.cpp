#include "velocityverlet.h"
#include "../box.h"


VelocityVerlet::VelocityVerlet(class Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{
    dt2 = dt / 2;
    ddt2 = dt * dt2;
}

void VelocityVerlet::next_step()
{
    /* Move particle one step
     */

    mat acc_old = box->accelerations;
    box->positions += box->velocities * dt + acc_old * ddt2;
    box->poteng = box->forcefield->eval_acc(box->positions, box->accelerations, box->potengs, true);
    box->velocities += (box->accelerations + acc_old) * dt2;
}
