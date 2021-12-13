#pragma once
#include "moves.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class Box*, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    bool reject_move;
    int n_in;
    double r_above, r_above2, v_in;
    class Particle* particle_out = nullptr;
};
