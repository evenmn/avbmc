#pragma once
#include <string>
#include <valarray>

#include "moves.h"


class Trans : public Moves
{
public:
    Trans(class Box*, double = 0.01);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    int i;
    double dx;
    std::valarray<double> pos_old;
};
