#pragma once
#include "avbmcin.h"
#include "avbmcout.h"


class AVBMC : public AVBMCIn, public AVBMCOut
{
public:
    AVBMC(class Box* box_in, const double r_below_in=0.95, const double r_above_in=3.0);
    void perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    bool move_in;
};
