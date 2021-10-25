#include "trans.h"
#include "../box.h"


Trans::Trans(class Box* box_in, double dx_in)
    : Moves(box_in) 
{
    dx = dx_in;
}

mat Trans::perform_move()
{
    /* Propose naive translation
     * move
     */

    // declare copy of positions
    mat pos_copy = box->positions;
    
    // pick particle to move
    int i = rng->next_int(box->npar);
    rowvec dr(3);
    for(double &i : dr){
        i = rng->next_double();
    }
    dr *= dx / norm(dr);
    pos_copy.row(i) += 2 * (rng->next_double() - 0.5) * dr; 
    return pos_copy;
}

double Trans::accept()
{
    /* Gives the ratio between the translation
     * probabilities Tij and Tji. For the 
     * naive case, Tij = Tji.
     */
    return 1.;
}
