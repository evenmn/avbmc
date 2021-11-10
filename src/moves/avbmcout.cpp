#include "avbmcout.h"
#include "../box.h"


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
    v_in = 4 * datum::pi * pow(r_above, 3)/3;
}

void AVBMCOut::perform_move(const int i)
{
    /* Remove a random particle from the bonded region
     * of particle i
     */

    
    if(box->npar < 2){
        box->sampler->du = 1e9;
    }
    else{
        // create local neighbor list of particle i
        std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
        n_in = neigh_listi.size();

        // pick particle to be removed
        int idx = box->rng->next_int(n_in);
        int j = neigh_listi[idx];

        box->sampler->du = box->potengs(j);
    }
}

double AVBMCOut::accept()
{
    /* 
     */
    return n_in * box->npar/ (v_in * (box->npar - 1));
}

void AVBMCOut::update_box(const int i)
{
    /*
     */
    box->npar --;
    box->forcefield->rm_distance_cross(i);
    box->positions.shed_row(i);
}
