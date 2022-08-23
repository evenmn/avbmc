/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <chrono>


class Constraint
{
public:
    Constraint(class Box *);
    Constraint(const Constraint &);
    virtual bool verify() = 0;
    virtual ~Constraint() = default;

    std::string label;
    double cum_time;
    unsigned int nreject;

protected:
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::duration<double> fsec;
    unsigned int cutoff_id, type1, type2;
    std::vector<std::vector<unsigned int> > neigh_list;
    class Box* box = nullptr;
};
