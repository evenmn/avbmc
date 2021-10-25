#pragma once
#include "moves.h"


class Trans : public Moves
{
public:
    Trans(class Box* box_in, double dx_in);
    mat perform_move();
    double accept();

private:
    double dx;
};
