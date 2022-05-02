#pragma once
#include <string>

#include "avbmcmolin.h"
#include "avbmcmolout.h"


class AVBMCMol : public AVBMCMolIn, public AVBMCMolOut
{
public:
    AVBMCMol(class System *, class Box *, std::vector<class Particle>, double = 0.95, double = 3.0, double = 1.3, bool = false, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_nsystemsize() override;
    std::string repr() override;

private:
    double r_above, r_below, r_inner;
    bool move_in, energy_bias, target_mol;
    class Box* box = nullptr;
};
