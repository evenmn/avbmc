/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <cmath>
#include <string>

#include "moves.h"


class AVBMCIn : virtual public Moves
{
public:
    AVBMCIn(class System *, class Box *, const std::string &, double = 0.9,
        double = 1.5, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool energy_bias;
    unsigned int n_in, particle_type;
    double v_in;
    std::string particle_label;
    class Box* box = nullptr;
};
