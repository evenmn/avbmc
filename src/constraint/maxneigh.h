/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <string>
#include <vector>

#include "constraint.h"


class MaxNeigh : public Constraint
{
public:
    MaxNeigh(class Box *, std::string, std::string, double, int);
    bool verify() override;

private:
    unsigned int nc;
};
