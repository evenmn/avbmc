#pragma once
#include <cmath>
#include "moves.h"


class TransMH : public Moves
{
public:
    TransMH(class Box* box_in, const double dx_in=0.01, const double Ddt_in=0.01);
    void perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    double dx, Ddt;
    rowvec eps, da;
};
