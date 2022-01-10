#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class AVBMCMolIn : virtual public Moves
{
public:
    AVBMCMolIn(class System *, class Box *, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool reject_move;
    int natom;
    double r_below, r_above, r_belowsq, r_abovesq, v_in, nmolavg;
    class Molecule* molecule_in = nullptr;
    class Box* box = nullptr;
};
