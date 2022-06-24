/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <valarray>

#include "boundary.h"


class Periodic : public Boundary
{
public:
    Periodic(class Box *, std::valarray<double>);
    void correct_position(unsigned int) override;
    void correct_distance(std::valarray<double> &) override;

private:
    unsigned int ndim;
    double volume;
    std::valarray<double> length;
};
