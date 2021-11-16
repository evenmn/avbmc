#include "trans.h"
#include "../box.h"


Trans::Trans(class Box* box_in, double dx_in)
    : Moves(box_in) 
{
    dx = dx_in;
}

void Trans::perform_move()
{
    /* Propose naive translation move of particle i
     */
    std::cout << "Translating particle" << std::endl;
    std::cout << box->npar << std::endl;

    // compute initial energy contribution from particle i
    i = box->rng->next_int(box->npar);
    double u0 = box->forcefield->comp_energy_par(box->particles, i);

    // move particle i
    std::valarray<double> dr(box->ndim);
    for(double &j : dr)
        j = rng->next_double();
    //dr /= std::sqrt(std::sum(std::pow(dr, 2)));  // normalize
    pos_old = box->particles[i]->r;
    box->particles[i]->r += 2 * (rng->next_double() - 0.5) * dx * dr;

    // compute new energy contribution from particle i
    double u1 = box->forcefield->comp_energy_par(box->particles, i);
    box->sampler->du = u1 - u0;
    box->poteng += box->sampler->du;

    //std::cout << "perform_move10" << std::endl;
    // compute new acceleration and potential energy
    //box->forcefield->tmp_to_state();
    //box->forcefield->update_distance_par(tmp_positions, i);
    //double u1 = box->forcefield->update_force_par(tmp_positions, i);
    //std::cout << "perform_move11" << std::endl;
    //box->sampler->du = u1 - box->poteng;
    //std::cout << "perform_move12" << std::endl;
    //box->boundary->update();
    //std::cout << "perform_move13" << std::endl;
}

double Trans::accept()
{
    /* Gives the ratio between the translation
     * probabilities Tij and Tji. For the 
     * naive case, Tij = Tji.
     */
    return 1.;
}

void Trans::reset()
{
    /* Set back to initial state if move is rejected
     */
    box->particles[i]->r = pos_old;
    box->poteng -= box->sampler->du;
}

/*
void Trans::update_box(const int i)
{
    cout << "update_box1" << endl;
    box->poteng += box->sampler->du;
    cout << "update_box2" << endl;
    box->positions = tmp_positions;
    cout << "update_box3" << endl;
    box->accelerations = box->forcefield->tmp_accelerations;
    cout << "update_box4" << endl;
    box->forcefield->state_to_tmp();
    cout << "update_box5" << endl;
    box->boundary->update();
    cout << "update_box6" << endl;
}
*/
