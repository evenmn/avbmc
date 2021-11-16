#pragma once
#include <cmath>
#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class Box* box_in, const double r_below_in=0.95, const double r_above_in=3.0);
    void perform_move();
    double accept();
    void reset();

private:
    int n_in;
    double r_below, r_above, r_above2, r_ratio, v_in;
    class Particle* particle_in;
};
