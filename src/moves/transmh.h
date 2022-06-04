/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <cmath>
#include <string>
#include <valarray>
#include "moves.h"


class TransMH : public Moves
{
public:
    TransMH(class System *, class Box *, double = 0.01, double = 0.01,
        const std::string & = "");
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool element_spec;
    unsigned int i, element_type;
    double dx, Ddt;
    std::valarray<double> eps, df;
    class Box* box = nullptr;
};
