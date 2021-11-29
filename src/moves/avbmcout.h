#pragma once
#include "moves.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class Box* box_in, const double r_above_in=3.0);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    bool not_accept;
    int n_in, npar_prev_prev;
    double r_above, r_above2, v_in;
    class Particle * particle_out;
};
