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
    unsigned int natom, neigh_id_above, neigh_id_below, neigh_id_inner;
    double r_below, r_above, r_inner, r_belowsq, r_abovesq, v_in, nmolavg, natom_inv;
    std::vector<int> npartype_old;
    std::vector<class Particle> particles, particles_old;
    class Box* box = nullptr;
};
