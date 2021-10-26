#pragma once
#include <cmath>
#include "moves.h"


class TransMH : public Moves
{
public:
    TransMH(class Box* box_in, double dx_in=0.01, double Ddt_in=0.01);
    mat perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    double dx, Ddt, u1;
    rowvec eps, a1;
};
