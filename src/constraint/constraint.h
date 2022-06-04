/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <string>
#include <vector>


class Constraint
{
public:
    Constraint(class Box *);
    Constraint(const Constraint &);
    virtual bool verify() = 0;
    virtual ~Constraint() = default;

    std::string label;

protected:
    unsigned int cutoff_id, type1, type2;
    std::vector<std::vector<unsigned int> > neigh_list;
    class Box* box = nullptr;
};
