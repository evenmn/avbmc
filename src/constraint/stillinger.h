/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <vector>
#include <string>
#include <valarray>

#include "constraint.h"


class Stillinger : public Constraint
{
public:
    Stillinger(class Box *, double = 2.0);
    Stillinger(const Stillinger &);  // copy constructor
    Stillinger& operator=(const Stillinger &other) // overloading equal operator
    {
        Stillinger tmp(other); // calling copy constructor
        swap(tmp);
        return *this;
    };
    void swap(Stillinger &other);  // swap two objects

    void set_criterion(std::string, std::string, double);
    void check_neigh_recu(int, std::valarray<char> &, std::valarray<char> &);
    bool verify() override;
    //double comp_volume();
    ~Stillinger();

private:
    unsigned int ntype; //, cutoff_id;
    double v_c, **r_csq_mat;
    std::vector<std::vector<int> > neigh_lists;
    std::valarray<char> in_cluster, checked;
};
