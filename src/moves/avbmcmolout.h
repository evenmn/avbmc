#pragma once
#include <string>
#include <vector>

#include "moves.h"


class AVBMCMolOut : virtual public Moves
{
public:
    AVBMCMolOut(class System *, class Box *, double = 3.0);
    AVBMCMolOut(class System *, class Box *, std::vector<Particle>, double = 3.0, double = 2.0);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    unsigned int natom;
    bool reject_move;
    double r_above, r_abovesq, v_in, nmolavg, r_max_inner;
    std::vector<class Particle> particles_old, molecule;
    class Box* box = nullptr;
};
