#pragma once
#include <cmath>
#include <string>
#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class Box* box_in, const double r_below_in=0.9, const double r_above_in=1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    std::vector<std::valarray<double> > rotate_molecule(std::vector<std::valarray<double> >);
    int n_in, type, npar_prev_prev;
    double r_below, r_above, r_below2, r_above2, v_in;
    std::string label;
    class Molecule* molecule_in = nullptr;
};
