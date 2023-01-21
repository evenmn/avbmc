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
    DistanceManager(class Box *, const double & = 1e-2, const std::size_t & = 100);
    unsigned int add_cutoff(const double &);
    unsigned int add_cutoff(const double &, const std::string &, const std::string &, const bool & = true);
    unsigned int add_cutoff(double **);
    void set(), reset();
    void set(const std::size_t &);
    void update_trans(const unsigned int &);
    void update_remove(const unsigned int &);
    void update_insert(const unsigned int &);

    std::vector<unsigned int> build_neigh_list(const std::vector<class Particle> &,
        const unsigned int &, const double &);
    std::vector<unsigned int> build_neigh_list(const unsigned int &, const double &);
    
    std::vector<unsigned int> mapid2vector;
    //std::vector<std::vector<std::vector<unsigned int> > > neigh_lists;
    std::vector<std::vector<std::vector<unsigned int> > > neigh_lists;
    std::vector<std::vector<double> > distance_mat;
    std::vector<std::vector<std::valarray<double> > > distance_cube;
    std::vector<double **> cutoff_mats;

private:
  void clear_neigh(const unsigned int &);
  void remove_neigh(const unsigned int &);
  void update_neigh(const unsigned int &, const unsigned int &, const double &);
  void update_neigh_k(const unsigned int &);
  void update_neigh_k(const unsigned int &, const unsigned int &, const unsigned int &, const double &);
  double normsq(const std::valarray<double> &);

  double cutoff_tol;
  unsigned int ncutoff, nmode, max_neigh;
  std::valarray<unsigned int> nmodes;
  std::vector<bool> mutuals;
  std::vector<double> cutoffs0, cutoffs1;
  std::vector<unsigned int> types1, types2, modes;

  std::vector<std::vector<std::vector<unsigned int> > > neigh_lists_old;
  std::vector<std::vector<double> > distance_mat_old;
  std::vector<std::vector<std::valarray<double> > > distance_cube_old;

  class Box *box = nullptr;
};


