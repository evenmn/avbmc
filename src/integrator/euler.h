/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include "integrator.h"


class Euler : public Integrator
{
public:
    Euler(class Box *, double = 0.001);
    double next_step() override;
};
