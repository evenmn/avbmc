#include "avbmcout.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"


AVBMCOut::AVBMCOut(Box* box_in, const double r_above_in)
    : Moves(box_in)
{
    /* Aggregate-volume-biased out move. 
     * This is a fictitious inter-box move, and works in the
     * grand canonical essemble only
     */

    // initialize r_above
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    not_accept = false;
}

void AVBMCOut::perform_move()
{
    /* Remove a random molecule from the bonded region
     * of another similar molecule
     */
    
    if(box->npar < 2){
        not_accept = true;
    }
    else{
        not_accept = false;
        // pick molecule to be removed
        int mol_idx = box->rng->choice(box->molecule_types->molecule_probs);
        molecule_out = box->molecule_types->construct_molecule(mol_idx);

        du = - box->forcefield->comp_energy_mol(box->particles, molecule_out);

        particles_old = box->particles;
        for(int atom : molecule_out->atoms_idx){
            box->particles.erase(box->particles.begin() + atom);
        }
        box->npar -= molecule_out->natom;

        /*
        // create local neighbor list of particle i

        int i = box->rng->next_int(box->npar);
        std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
        n_in = neigh_listi.size();

        if(n_in > 0){
            // pick particle to be removed
            int neigh_idx = box->rng->next_int(n_in);
            int j = neigh_listi[neigh_idx];  // particle to be removed
            du = - box->forcefield->comp_energy_par(box->particles, j);
            box->poteng += du;
            particle_out = box->particles[j];
            box->particles.erase(box->particles.begin() + j);
            box->npar --;
        }
        else{
            not_accept = true;
        }
        */
    }
}

double AVBMCOut::accept(double temp, double chempot)
{
    /* Return the acceptance probability of move
     */
    if(not_accept){
        return 0;
    }
    else{
        double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar+1);
        return n_in * box->npar / (v_in * (box->npar - 1)) * std::exp(-(du+chempot+dw)/temp);
    }
}

void AVBMCOut::reset()
{
    /* Set back to old state if
     * move is rejected
     */
    box->npar += molecule_out->natom;
    box->poteng -= du;
    box->particles = particles_old;
}
