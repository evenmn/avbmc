#pragma once
#include <string>

#include "avbmcmolin.h"
#include "avbmcmolout.h"


class AVBMCMol : public AVBMCMolIn, public AVBMCMolOut
{
public:
    AVBMCMol(class System *, class Box *, std::vector<class Particle>, double, double = 0.95, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    double r_above, r_below;
    bool move_in;
    class Box* box = nullptr;
};
