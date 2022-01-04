#pragma once
#include <string>
#include <memory>

#include "moves.h"


class AVBMCOut : virtual public Moves
{
public:
    //AVBMCOut(class System *, class Box *, double = 3.0);
    AVBMCOut(class System *, std::shared_ptr<class Box>, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    bool reject_move;
    int n_in;
    double r_above, r_abovesq, v_in;
    std::shared_ptr<class Particle> particle_out = nullptr;
    //class Box* box = nullptr;
    std::shared_ptr<class Box> box = nullptr;
};
