#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <valarray>

#include "moves.h"


class AVBMCMolInRes : virtual public Moves
{
public:
    AVBMCMolInRes(class System *, class Box *, std::vector<class Particle>,
    double, double = 0.9, double = 1.5, bool = true, bool = false);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool reject_move, energy_bias, target_mol;
    unsigned int natom;
    double r_below, r_above, r_max_inner, r_belowsq, r_abovesq, v_in, nmolavg, natom_inv;
    std::vector<class Particle> particles;
    std::vector<class Particle> particles_old;
    class Box* box = nullptr;
};
