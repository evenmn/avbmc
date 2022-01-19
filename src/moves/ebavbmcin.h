#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class EBAVBMCIn : virtual public Moves
{
public:
    EBAVBMCIn(class System *, class Box *, std::string, double = 0.9, double = 1.5);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    unsigned int n_in, particle_type;
    double r_below, r_above, r_belowsq, r_abovesq, v_in;
    std::string particle_label;
    class Box* box = nullptr;
};
