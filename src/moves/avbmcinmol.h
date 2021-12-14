#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class AVBMCInMol : virtual public Moves
{
public:
    AVBMCInMol(class Box*, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    bool reject_move;
    int n_in, natom;
    std::string label;
    double r_below, r_above, r_below2, r_above2, v_in;
    class Molecule* molecule_in = nullptr;
};
