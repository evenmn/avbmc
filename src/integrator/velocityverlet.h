/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include "integrator.h"


class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet(class Box *, double = 0.01);
    double next_step() override;

private:
    double dt2;
    double ddt2;
};
