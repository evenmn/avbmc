#pragma once
#include "moves.h"


class AVBMCOut : public Moves
{
public:
    AVBMCOut(class Box* box_in, const double p_bias_in, const double r_above_in);
    void perform_move(const int i);
    double accept();
    void update_box(const int i);

private:
    double p_bias, r_above, v_in, p_ratio;
};
