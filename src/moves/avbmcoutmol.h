#pragma once
#include "moves.h"


class AVBMCOutMol : virtual public Moves
{
public:
    AVBMCOutMol(class Box*, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    bool not_accept;
    int n_in;
    double r_above, r_above2, v_in;
    std::vector<class Particle *> particles_old;
    class Molecule* molecule_out = nullptr;
};
