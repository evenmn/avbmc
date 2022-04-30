#pragma once
#include <string>

#include "avbmcin.h"
#include "avbmcout.h"


class AVBMC : public AVBMCIn, public AVBMCOut
{
public:
    AVBMC(class System *, class Box *, std::string, double = 0.95, double = 3.0);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_nsystemsize() override;
    std::string repr() override;

private:
    bool move_in;
    double r_below, r_above;
    class Box* box = nullptr;
    std::string particle_label;
};
