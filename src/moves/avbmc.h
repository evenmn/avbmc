#pragma once
#include <string>

#include "avbmcin.h"
#include "avbmcout.h"


class AVBMC : public AVBMCIn, public AVBMCOut
{
public:
    AVBMC(class System *, class Box *, std::string, double = 0.95, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    bool move_in;
    double r_below, r_above;
    class Box* box = nullptr;
};
