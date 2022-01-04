#pragma once
#include <cmath>
#include <string>
#include <memory>

#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    //AVBMCIn(class System *, class Box *, double = 0.9, double = 1.5);
    AVBMCIn(class System *, std::shared_ptr<class Box>, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    //std::vector<std::valarray<double> > rotate_molecule(std::vector<std::valarray<double> >);
    int n_in, type;
    double r_below, r_above, r_belowsq, r_abovesq, v_in;
    std::string label_in;
    //class Box* box = nullptr;
    std::shared_ptr<class Box> box = nullptr;
};
