#pragma once
#include "moves.h"


class Trans : public Moves
{
public:
    Trans(class Box* box_in, double dx_in=0.01);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    int i;
    double dx;
    std::valarray<double> pos_old;
};
