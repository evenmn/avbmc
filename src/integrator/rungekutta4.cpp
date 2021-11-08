#include "rungekutta4.h"
#include "../box.h"


RungeKutta4::RungeKutta4(class Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}

void RungeKutta4::next_step()
{
    /* Move particle one step
     */

    // store old positions and velocities
    pos_old = box->positions;
    vel_old = box->velocities;

    // calculate K1
    K1x = box->velocities * dt;
    K1v = box->accelerations * dt;

    // calculate K2
    pos_new = pos_old + K1x / 2;
    vel_new = vel_old + K1v / 2;
    box->forcefield->eval_acc(pos_new, acc_new, potengs, false);
    K2x = vel_new * dt;
    K2v = acc_new * dt;

    // calculate K3
    pos_new = pos_old + K2x / 2;
    vel_new = vel_old + K2v / 2;
    box->forcefield->eval_acc(pos_new, acc_new, potengs, false);
    K3x = vel_new * dt;
    K3v = acc_new * dt;

    // calculate K4
    pos_new = pos_old + K3x;
    vel_new = vel_old + K3v;
    box->forcefield->eval_acc(pos_new, acc_new, potengs, false);
    K4x = vel_new * dt;
    K4v = acc_new * dt;
    
    // move
    box->positions = pos_old + (K1x + 2 * (K2x + K3x) + K4x) / 6;
    box->velocities = vel_old + (K1v + 2 * (K2v + K3v) + K4v) / 6;
    box->poteng = box->forcefield->eval_acc(box->positions, box->accelerations, box->potengs, true);
}
