#pragma once
//#include <cmath>
#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class Box* box_in, const double r_below_in=0.9, const double r_above_in=1.5);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    int n_in, type;
    double r_below, r_above, r_below2, r_above2, v_in;
    std::string label;
    //class Particle* particle_in = nullptr;
};
