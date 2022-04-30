#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class System *, class Box *, std::string, double = 0.9, double = 1.5);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_nsystemsize() override;
    std::string repr() override;

private:
    unsigned int n_in, particle_type;
    double r_below, r_above, r_belowsq, r_abovesq, v_in;
    std::string particle_label;
    class Box* box = nullptr;
};
