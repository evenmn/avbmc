/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>

#include "moves.h"
#include "../particle.h"


class AVBMCOut : virtual public Moves
{
public:
    AVBMCOut(class System *, class Box *, const std::string &, double = 3.0,
        bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool move_performed, energy_bias;
    unsigned int n_in, particle_type;
    std::string particle_label;
    double v_in;
    Particle particle_out = Particle("", {0});
    class Box* box = nullptr;
};
