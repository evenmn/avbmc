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
    v_in = 4 * datum::pi * pow(r_above, 3)/3;

    type = 0;
    chem_symbol = "Ar";
    mass = 1.;
}

void AVBMCIn::perform_move(const int i)
{
    /* Remove a random particle from the bonded region
     * of particle i
     */

    // create local neighbor list of particle i
    std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
    n_in = neigh_listi.size();

    // propose new particle position
    rowvec posi = box->positions.row(i);
    rowvec dr(box->ndim, fill::zeros);
    while(norm(dr) > r_ratio){
        for(int d=0; d<box->ndim; d++){
            dr(d) = box->rng->next_double();
        }
    }
    posj = posi + r_above * dr;

    // compute du
    box->sampler->du = box->forcefield->comp_force_par(posj, accj);
}

double AVBMCIn::accept()
{
    /* 
     */
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1));
}

void AVBMCIn::update_box(const int i)
{
    /*
     */
    box->npar ++;
    box->positions.insert_rows(box->npar-2, posj);
    box->forcefield->add_distance_cross(box->positions);
    box->potengs.insert_cols(box->npar-2, box->sampler->du);
    box->particle_types.push_back(type);
    box->chem_symbols.push_back(chem_symbol);
    box->particle_masses.push_back(mass);
}
