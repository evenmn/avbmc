#pragma once
#include <string>

#include "moves.h"
#include "../particle.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class System *, class Box *, const std::string &, double = 3.0);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_nsystemsize() override;
    std::string repr() override;

private:
    bool reject_move;
    unsigned int n_in;
    std::string particle_label;
    double r_above, r_abovesq, v_in;
    Particle particle_out = Particle("", {0});
    class Box* box = nullptr;
};
