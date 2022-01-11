#pragma once
#include <string>

#include "moves.h"
#include "../particle.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class System *, class Box *, std::string, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool reject_move;
    unsigned int n_in;
    std::string particle_label;
    double r_above, r_abovesq, v_in;
    Particle particle_out = Particle("", {0});
    class Box* box = nullptr;
};
