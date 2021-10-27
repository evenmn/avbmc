#pragma once
#include "moves.h"


class Trans : public Moves
{
public:
    Trans(class Box* box_in, double dx_in=0.01);
    void perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    double dx, u1;
    rowvec eps, a1;
};
