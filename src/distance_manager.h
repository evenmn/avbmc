/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <valarray>


class DistanceManager
{
public:
    DistanceManager(class Box *, double = 1e-2);
    unsigned int add_cutoff(double);
    unsigned int add_cutoff(double, std::string, std::string, bool = true);
    unsigned int add_cutoff(double **);
    void set(), reset();
    void update_trans(unsigned int);
    void update_remove(unsigned int);
    void update_insert(unsigned int);
    
    std::vector<std::vector<std::vector<unsigned int> > > neigh_lists;
    std::vector<std::vector<double> > distance_mat;
    std::vector<std::vector<std::valarray<double> > > distance_cube;
    std::vector<double **> cutoff_mats;

private:
  void clear_neigh(unsigned int);
  void remove_neigh(unsigned int);
  void update_neigh(unsigned int, unsigned int, double);
  double normsq(std::valarray<double>);

  double cutoff_tol;
  unsigned int ncutoff;
  std::vector<bool> mutuals;
  std::vector<double> cutoffs;
  std::vector<unsigned int> types1, types2, modes;

  std::vector<std::vector<std::vector<unsigned int>>> neigh_lists_old;
  std::vector<std::vector<double>> distance_mat_old;
  std::vector<std::vector<std::valarray<double>>> distance_cube_old;

  class Box *box = nullptr;
};


