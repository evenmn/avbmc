#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <valarray>

#include "moves.h"


class AVBMCMolIn : virtual public Moves
{
public:
    AVBMCMolIn(class System *, class Box *, double = 0.9, double = 1.5);
    AVBMCMolIn(class System *, class Box *, std::vector<std::string>, std::vector<std::valarray<double> >, double, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool reject_move;
    unsigned int natom;
    double r_below, r_above, r_max_inner, r_belowsq, r_abovesq, v_in, nmolavg;
    std::vector<std::string> molecule_elements;
    std::vector<std::valarray<double> > atom_positions;
    std::vector<int> atom_types;
    class Box* box = nullptr;
};
