#include "avbmcout.h"
#include "../box.h"
#include "../particle.h"


AVBMCOut::AVBMCOut(class Box* box_in, const double r_above_in)
    : Moves(box_in)
{
    /* Aggregate-volume-biased out move. 
     * This is a fictitious inter-box move, and work is the
     * grand canonical essemble only
     */

    // initialize r_above
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    v_in = 4 * pi * pow(r_above, 3)/3;
    not_accept = false;
}

void AVBMCOut::perform_move()
{
    /* Remove a random particle from the bonded region
     * of particle i
     */
    
    if(box->npar < 2){
        not_accept = true;
    }
    else{
        not_accept = false;
        // create local neighbor list of particle i
        int i = box->rng->next_int(box->npar);
        std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
        n_in = neigh_listi.size();

        if(n_in > 0){
            // pick particle to be removed
            int neigh_idx = box->rng->next_int(n_in);
            int j = neigh_listi[neigh_idx];  // particle to be removed
            box->sampler->du = - box->forcefield->comp_energy_par(box->particles, j);
            particle_out = box->particles[j];
            box->particles.erase(box->particles.begin() + j);
            box->npar --;
        }
        else{
            not_accept = true;
        }
    }
}

double AVBMCOut::accept(double temp, double chempot)
{
    /* Return the pre-exponential factor before
     * Boltzmann ratio
     */
    if(not_accept){
        return 0;
    }
    else{
        return n_in * box->npar / (v_in * (box->npar - 1)) * std::exp(-chempot/temp);
    }
}

void AVBMCOut::reset()
{
    box->npar ++;
    box->poteng -= box->sampler->du;
    box->particles.push_back(particle_out);
}
