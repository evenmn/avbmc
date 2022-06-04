/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>


class Integrator
{
public:
    Integrator(class Box *, double);
    virtual double next_step() = 0;
    virtual ~Integrator() = default;

    double dt;

protected:
    class Box* box = nullptr;
};
