#pragma once
#include <string>
#include <vector>

#include "moves.h"


class AVBMCMolOutRes : virtual public Moves
{
public:
    AVBMCMolOutRes(class System *, class Box *, std::vector<Particle>,
                   double = 3.0, double = 2.0, bool = true, bool = false);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    unsigned int natom;
    bool reject_move, energy_bias, target_mol;
    double r_above, r_abovesq, v_in, nmolavg, r_max_inner, natom_inv;
    std::vector<class Particle> particles_old, molecule;
    class Box* box = nullptr;
};
