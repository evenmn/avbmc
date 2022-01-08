#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class System *, class Box *, std::string, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    int n_in, type;
    double r_below, r_above, r_belowsq, r_abovesq, v_in;
    std::string label_in;
    class Box* box = nullptr;
};
