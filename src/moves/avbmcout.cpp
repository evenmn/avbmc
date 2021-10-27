#include "avbmcout.h"
#include "../box.h"


AVBMCOut::AVBMCOut(class Box* box_in, const double p_bias_in, const double r_above_in)
    : Moves(box_in)
{
    // initialize p_bias
    p_bias = p_bias_in;
    p_ratio = (1-p_bias)/p_bias;

    // initialize r_above
    r_above = r_above_in;
    v_in = 4 * datum::pi * pow(r_above, 3)/3;
}

void AVBMCOut::perform_move(const int i)
{
    /* Remove particle i
     */
    mat pos_copy = box->positions;
    pos_copy.shed_row(i);

    box->sampler->du = box->potengs(i);
    //sampler->da = box->accelerations.row(i);
}

double AVBMCOut::accept()
{
    /* 
     */
    double v_out = box->boundary->comp_volume() - v_in;
    return p_ratio * v_in / v_out;
}

void AVBMCOut::update_box(const int i)
{
    /*
     */
    box->forcefield->rm_distance_cross(i);
    box->positions.shed_row(i);
}
