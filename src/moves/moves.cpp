#include "moves.h"
#include "../box.h"

Moves::Moves(class Box* box_in)
{
    box = box_in;
    rng = box->rng;
}

bool Moves::check_stillinger()
{
    /* Check Stillinger criterion
     */
    return false;
}

double Moves::get_cluster_volume()
{
    /* Get volume of cluster
     */
    return 0.;
}

//Moves::~Moves() {}
