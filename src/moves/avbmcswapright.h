/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>

#include "moves.h"
#include "avbmcin.h"
#include "avbmcout.h"


class AVBMCSwapRight : public AVBMCIn, public AVBMCOut
{
public:
    AVBMCSwapRight(class System *, class Box *, class Box *,
        const std::string &, double = 0.95, double = 3.0, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool energy_bias;
    double r_below, r_above;
    class Particle* particle_out = nullptr;
    class Box* box1 = nullptr;
    class Box* box2 = nullptr;
};
