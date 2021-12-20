#pragma once
#include <string>

#include "moves.h"


class AVBMCSwapRight : virtual public Moves
{
public:
    AVBMCOut(class System *, class Box *, class Box *, double = 0.95, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    bool reject_move;
    int n_in;
    double r_above, r_abovesq, v_in;
    class Particle* particle_out = nullptr;
    class Box* box1 = nullptr;
    class Box* box2 = nullptr;
};
