/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>

#include "moves.h"
#include "avbmcmolin.h"
#include "avbmcmolout.h"


class AVBMCMolSwapRight : public AVBMCMolIn, public AVBMCMolOut
{
public:
    AVBMCMolSwapRight(class System *, class Box *, class Box *,
        std::vector<Particle>, double = 0.95, double = 3.0, double = 1.5,
        bool = false, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool energy_bias, target_mol;
    double r_below, r_above, r_inner;
    class Particle* particle_out = nullptr;
    class Box* box1 = nullptr;
    class Box* box2 = nullptr;
};
