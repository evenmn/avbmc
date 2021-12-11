#pragma once
#include <cmath>
#include "moves.h"


class AVBMCInMol : virtual public Moves
{
public:
    AVBMCInMol(class Box*, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    std::vector<std::valarray<double> > rotate_molecule(std::vector<std::valarray<double> >);
    int n_in, natom;
    double r_below, r_above, r_below2, r_above2, v_in;
    std::string label;
    class Molecule* molecule_in = nullptr;
};
