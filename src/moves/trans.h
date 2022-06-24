/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>
#include <valarray>

#include "moves.h"


class Trans : public Moves
{
public:
    Trans(class System *, class Box *, double = 0.01, const std::string & = "");
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    bool element_spec;
    unsigned int i, element_type;
    double dx;
    std::valarray<double> pos_old;
    class Box* box = nullptr;
};
