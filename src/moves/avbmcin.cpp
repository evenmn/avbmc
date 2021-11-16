#include "avbmcin.h"
#include "../box.h"


AVBMCIn::AVBMCIn(class Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(box_in)
{
    /* Aggregate-volume-biased out move. 
     * This is a fictitious inter-box move, and work is the
     * grand canonical essemble only
     */

    r_below = r_below_in;
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    r_ratio = r_below / r_above;
    v_in = 4 * pi * pow(r_above, 3)/3;

    particle_in->type = 0;
    particle_in->mass = 1.;
    particle_in->label = "Ar";
}

void AVBMCIn::perform_move()
{
    /* Remove a random particle from the bonded region
     * of particle i
     */

    std::cout << "Inserting particle" << std::endl;
    std::cout << box->npar << std::endl;

    // create local neighbor list of particle i
    int i = box->rng->next_int(box->npar);
    std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
    n_in = neigh_listi.size();

    auto norm = [] (std::valarray<double> x) -> double { return std::sqrt(std::pow(x, 2).sum()); };

    // construct new particle
    std::valarray<double> dr(box->ndim);
    while(norm(dr) < r_ratio){
        for(int d=0; d<box->ndim; d++){
            dr[d] = box->rng->next_double();
        }
    }
    particle_in->r = box->particles[i]->r + r_above * dr;
    box->particles.push_back(particle_in);
    box->npar ++;

    // compute du
    box->sampler->du = box->forcefield->comp_energy_par(box->particles, box->npar - 1);
    box->poteng += box->sampler->du;

    /*
    cout << "perform_move1" << endl;
    tmp_positions = box->positions;
    cout << "perform_move2" << endl;
    posj = tmp_positions.row(i) + r_above * dr;
    cout << "perform_move3" << endl;
    tmp_positions.insert_rows(box->npar - 1, posj);
    cout << "perform_move4" << endl;

    // compute du
    box->forcefield->add_distance_par(tmp_positions);
    cout << "perform_move5" << endl;
    double u1 = box->forcefield->update_force_par(tmp_positions, box->npar);
    cout << "perform_move6" << endl;
    box->sampler->du = u1 - box->poteng;
    cout << "perform_move7" << endl;
    box->boundary->update();
    cout << "perform_move8" << endl;
    */
}

double AVBMCIn::accept()
{
    /* 
     */
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1));
}

void AVBMCIn::reset()
{
    box->npar --;
    box->poteng -= box->sampler->du;
    box->particles.erase(box->particles.begin() + box->npar);
}

/*
void AVBMCIn::update_box(const int i)
{
    box->npar ++;
    box->poteng += box->sampler->du;
    std::cout << "npar: " << box->npar << std::endl;
    box->positions = tmp_positions;
    std::cout << "par rows: " << box->positions.n_rows << std::endl;
    box->accelerations = box->forcefield->tmp_accelerations;
    std::cout << "acc rows: " << box->accelerations.n_rows << std::endl;
    box->forcefield->state_to_tmp();
    box->boundary->update();
}
*/
