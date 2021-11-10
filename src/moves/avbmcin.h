#pragma once
#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class Box* box_in, const double r_below_in=0.95, const double r_above_in=3.0);
    void perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    std::string chem_symbol;
    int n_in, type;
    double r_below, r_above, r_above2, r_ratio, v_in, mass;
    rowvec posj, accj;
};
