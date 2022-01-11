#pragma once
#include <string>
#include <vector>

#include "moves.h"


class AVBMCMolOut : virtual public Moves
{
public:
    AVBMCMolOut(class System *, class Box *, double = 3.0);
    AVBMCMolOut(class System *, class Box *, std::vector<std::string>, double = 3.0, double = 2.0);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool reject_move;
    double r_above, r_abovesq, v_in, nmolavg;
    std::vector<class Particle> particles_old;
    class Molecule* molecule_out = nullptr;
    class Box* box = nullptr;
};
