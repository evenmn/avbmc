#pragma once
#include <string>
#include <vector>

#include "moves.h"


class AVBMCOutMol : virtual public Moves
{
public:
    AVBMCOutMol(class System *, class Box *, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    bool reject_move;
    double r_above, r_abovesq, v_in, nmolavg;
    std::vector<class Particle *> particles_old;
    class Molecule* molecule_out = nullptr;
    class Box* box = nullptr;
};
