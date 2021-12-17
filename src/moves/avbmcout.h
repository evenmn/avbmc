#pragma once
#include <string>

#include "moves.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class Box*, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    bool reject_move;
    int n_in;
    double r_above, r_abovesq, v_in;
    class Particle* particle_out = nullptr;
};
