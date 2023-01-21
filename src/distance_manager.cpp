/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-23)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   The distance manager is used to manage the distance matrix and the neighbor
   lists. This includes storing the distance matrix, the relative coordinates
   between particles and all neighbor lists, and updating the distance matrix
   and neighbor lists when a particle is added, moved or removed. All parts
   of the code relying on the distance between particles (forcefield, moves
   etc..) call the distance manager - distance should never be computed outside
   this class!
---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <valarray>
#include <algorithm>  // for removing vector element by value

#include "distance_manager.h"
#include "box.h"
#include "forcefield/forcefield.h"
#include "boundary/boundary.h"


/* ----------------------------------------------------------------------------
   Distance manager class constructor. 'cutoff_tol_in' is the tolerance when
   checking if two cutoffs are equivalent and should share neighbor lists.
   'max_neigh_in' is the max number of neighbors an atom can take
---------------------------------------------------------------------------- */

DistanceManager::DistanceManager(Box* box_in, const double &cutoff_tol_in,
    const std::size_t &max_neigh_in)
{
    ncutoff = 0;
    nmode = 3;
    box = box_in;
    cutoff_tol = cutoff_tol_in;
    max_neigh = max_neigh_in;
    nmodes.resize(nmode, 0);
}


/* ----------------------------------------------------------------------------
   Add a cutoff to the list of cutoffs in order to make the distance manager
   manage neighbor lists. This cutoff applies to all components. Returning the
   cutoff-ID, which has to be used when extracting the correct neighbor list.
---------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(const double &rc)
{
    double rcsq;
    unsigned int k, mode, neigh_id;

    neigh_id = ncutoff;
    mode = 0;
    rcsq = rc * rc;

    // check if a similar cutoff exists
    for (k=0; k<ncutoff; k++) {
      if (modes[k] == mode && fabs(cutoffs0[mapid2vector[k]] - rcsq) < cutoff_tol) {
        return k;
      }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    modes.push_back(mode);
    mapid2vector.push_back(nmodes[mode]); 
    nmodes[mode] ++;
    cutoffs0.push_back(rcsq);

    // all particles
    neigh_lists.push_back({});
    neigh_lists[neigh_id].resize(box->npar);
    update_neigh_k(neigh_id);
    return neigh_id;
}


/* ----------------------------------------------------------------------------
   Add a cutoff 'rc' to the list of cutoffs in order to make the distance
   manager manage neighbor lists. Here, the cutoff is restricted for two
   components 'label1' and 'label2'. If 'mutual' is true, both label1 and 
   label2's neighbor lists will be updated. Returning the cutoff-ID, which has
   to be used when extracting the correct neighbor list.
---------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(const double &rc, const std::string &label1,
                                         const std::string &label2, const bool &mutual)
{
    double rcsq;
    unsigned int k, vecid, type1, type2, mode, neigh_id;

    neigh_id = ncutoff;
    mode = 1;
    rcsq = rc * rc;
    type1 = box->forcefield->label2type.at(label1);
    type2 = box->forcefield->label2type.at(label2);

    // check if a similar cutoff exists
    for (k=0; k<ncutoff; k++) {
        if (modes[k] == mode) {
        vecid = mapid2vector[k];
            if (
                types1[vecid] == type1 && 
                types2[vecid] == type2 &&
                mutuals[vecid] == mutual &&
                fabs(cutoffs1[vecid] - rcsq) < cutoff_tol
            ) {
                return k;
            }
        }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    modes.push_back(mode);
    mapid2vector.push_back(nmodes[mode]); 
    nmodes[mode] ++;
    types1.push_back(type1);
    types2.push_back(type2);
    mutuals.push_back(mutual);
    cutoffs1.push_back(rcsq);
    neigh_lists.push_back({});
    neigh_lists[neigh_id].resize(box->npar);
    update_neigh_k(neigh_id);
    return neigh_id;
}


/* ----------------------------------------------------------------------------
   Add a cutoff matrix 'rc' of shape (ntype, ntype) and create a combined
   neighbor list. Returning the cutoff-ID, which has to be used when extracting
   the correct neighbor list.
---------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double **rc)
{
    unsigned int mode, neigh_id;

    neigh_id = ncutoff;
    mode = 2;
    ncutoff ++;
    modes.push_back(mode);
    mapid2vector.push_back(nmodes[mode]); 
    nmodes[mode] ++;
    cutoff_mats.push_back(rc);
    neigh_lists.push_back({});
    neigh_lists[neigh_id].resize(box->npar);
    update_neigh_k(neigh_id);
    return neigh_id;
}


/* ----------------------------------------------------------------------------
   Clear neighbor lists of particle 'i' and clear it from the other neighbor
   lists. This is done when a particle is moved.

   TODO: This does not necessarily have to be done: perform checks first;
---------------------------------------------------------------------------- */

void DistanceManager::clear_neigh(const unsigned int &i) 
{
  std::vector<unsigned int>::iterator position;

  for (unsigned int j = 0; j < ncutoff; j++) {
    neigh_lists[j][i].clear(); // clear neighbor list of particle i
    for (unsigned int k = 0; k < neigh_lists[j].size(); k++) {
      position =
          std::find(neigh_lists[j][k].begin(), neigh_lists[j][k].end(), i);
      if (position != neigh_lists[j][k].end()) {
        neigh_lists[j][k].erase(position);
      }
    }
  }
}


/* ----------------------------------------------------------------------------
   Remove neighbor lists of particle 'i' and clear it from the other neighbor
   lists. This has to be done when a particle is removed.
---------------------------------------------------------------------------- */

void DistanceManager::remove_neigh(const unsigned int &i) {
  std::vector<unsigned int>::iterator position;

  for (unsigned int j = 0; j < ncutoff; j++) {
    neigh_lists[j].erase(neigh_lists[j].begin() + i);
    for (unsigned int k = 0; k < neigh_lists[j].size(); k++) {
      position =
          std::find(neigh_lists[j][k].begin(), neigh_lists[j][k].end(), i);
      if (position != neigh_lists[j][k].end()) {
        neigh_lists[j][k].erase(position);
      }
      for (unsigned int l = 0; l < neigh_lists[j][k].size(); l++) {
        if (neigh_lists[j][k][l] > i) {
          neigh_lists[j][k][l] -= 1;
        }
      }
    }
  }
}


/* ----------------------------------------------------------------------------
   Update neighbor list 'k' of a particle pair 'i' and 'j' with respect to the
   distance (squared), 'rij'.
---------------------------------------------------------------------------- */

void DistanceManager::update_neigh_k(const unsigned int &i, const unsigned int &j,
    const unsigned int &k, const double &rij)
{
    std::size_t typei = box->particles[i].type;
    std::size_t typej = box->particles[j].type;
    std::size_t vecid = mapid2vector[k];

    if (modes[k] == 0) {
        if (rij < cutoffs0[vecid]) {
            neigh_lists[k][i].push_back(j);
            neigh_lists[k][j].push_back(i);
        }
    }
    else if (modes[k] == 1) {
        if (types1[vecid]==typei && types2[vecid]==typej && rij < cutoffs1[vecid]) {
            neigh_lists[k][i].push_back(j);
            if (mutuals[vecid]) {
                neigh_lists[k][j].push_back(i);
            }
        }
    }
    else {
        if (rij < cutoff_mats[vecid][typei][typej]) {
            neigh_lists[k][i].push_back(j);
            neigh_lists[k][j].push_back(i);
        }
    }
}


/* ----------------------------------------------------------------------------
   Update all elements of neighbor list k
---------------------------------------------------------------------------- */

void DistanceManager::update_neigh_k(const unsigned int &k)
{
    for (std::size_t i=0; i<box->npar; i++) {
        for (std::size_t j=0; j<i; j++) {
            update_neigh_k(i, j, k, distance_mat[i][j]);
        }
    }
}


/* ----------------------------------------------------------------------------
   Update neighbor lists of a particle pair 'i' and 'j' with respect to the
   distance (squared), 'rij'. This is done when a particle is moved or inserted
---------------------------------------------------------------------------- */

void DistanceManager::update_neigh(const unsigned int &i, const unsigned int &j, const double &rij)
{
    /*
    unsigned int k, l, m, typei, typej;
    typei = box->particles[i].type;
    typej = box->particles[j].type;
    l = m = 0;
    for (k=0; k<ncutoff; k++) {
        if (modes[k] == 0) {
            if (rij < cutoffs[m]) {
                neigh_lists[k][i].push_back(j);
                neigh_lists[k][j].push_back(i);
            }
            m++;
        }
        else if (modes[k] == 1) {
            if (types1[l]==typei && types2[l]==typej && rij < cutoffs[m]) {
                neigh_lists[k][i].push_back(j);
                if (mutuals[l]) {
                    neigh_lists[k][j].push_back(i);
                }
            }
            l++;
            m++;
        }
        else {
            if (rij < cutoff_mats[k-m][typei][typej]) {
                neigh_lists[k][i].push_back(j);
                neigh_lists[k][j].push_back(i);
            }
        }
    }
    */
    for (unsigned int k=0; k<ncutoff; k++) {
        update_neigh_k(i, j, k, rij);
    }
}


/* ----------------------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
---------------------------------------------------------------------------- */

double DistanceManager::normsq(const std::valarray<double> &array)
{
    double norm = 0.;
    for (const double &element : array) {
        norm += element * element;
    }
    return norm;
}


/* ----------------------------------------------------------------------------
   ff
---------------------------------------------------------------------------- */

std::vector<unsigned int> DistanceManager::build_neigh_list(const unsigned int &i, const double &rsq)
{
    //double rijsq;
    //std::valarray<double> ri = box->particles[i].r;
    std::vector<unsigned int> neigh_list;

    for(unsigned int j=0; j<i; j++){
        //rijsq = normsq(box->particles[j].r - ri);
        if(distance_mat[i][j] < rsq){
            neigh_list.push_back(j);
        }
    }
    for(unsigned int j=i+1; j<box->npar; j++){
        //rijsq = normsq(box->particles[j].r - ri);
        if(distance_mat[i][j] < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ----------------------------------------------------------------------------
   ff
---------------------------------------------------------------------------- */

std::vector<unsigned int> DistanceManager::build_neigh_list(
    const std::vector<Particle> &particles, const unsigned int &i, const double &rsq)
{
    std::size_t npar = particles.size();
    std::valarray<double> ri = particles[i].r;
    std::vector<unsigned int> neigh_list;
    for (std::size_t j=0; j<i; j++) {
        double rijsq = std::pow(particles[j].r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for (std::size_t j=i+1; j<npar; j++) {
        double rijsq = std::pow(particles[j].r - ri, 2).sum();
        if (rijsq < rsq) {
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ----------------------------------------------------------------------------
   Store old matrices. Here the vectors have to be deep copied
---------------------------------------------------------------------------- */

void DistanceManager::set()
{
    distance_mat_old = distance_mat;
    distance_cube_old = distance_cube;
    neigh_lists_old = neigh_lists;
}

void DistanceManager::set(const std::size_t &i)
{
    for (std::size_t j=0; j<box->npar; j++) {
        distance_mat_old[i][j] = distance_mat[i][j];
        distance_mat_old[j][i] = distance_mat[j][i];
        distance_cube_old[i][j] = distance_cube[i][j];
        distance_cube_old[j][i] = distance_cube[j][i];
    }
    neigh_lists_old = neigh_lists;
}


/* ----------------------------------------------------------------------------
   Reset new matrices. Only need a shallow copy here
---------------------------------------------------------------------------- */

void DistanceManager::reset()
{
    distance_mat = distance_mat_old;
    distance_cube = distance_cube_old;
    neigh_lists = neigh_lists_old;
}


/* ----------------------------------------------------------------------------
   Update distance matrix due to translational move of a particle 'i'. This 
   corresponds to updating the i'th row and column in the distance matrix.
---------------------------------------------------------------------------- */

void DistanceManager::update_trans(const unsigned int &i)
{
    double rij;
    std::valarray<double> delij;

    clear_neigh(i);
    for (std::size_t j=0; j<i; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        box->boundary->correct_distance(delij);
        rij = normsq(delij);
        distance_mat[i][j] = rij;
        distance_mat[j][i] = rij;
        distance_cube[i][j] = delij;
        distance_cube[j][i] = -delij;
        update_neigh(i, j, rij);
    }
    for (std::size_t j=i+1; j<box->npar; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        box->boundary->correct_distance(delij);
        rij = normsq(delij);
        distance_mat[i][j] = rij;
        distance_mat[j][i] = rij;
        distance_cube[i][j] = delij;
        distance_cube[j][i] = -delij;
        update_neigh(i, j, rij);
    }
}


/* ----------------------------------------------------------------------------
   Update distance matrix after a deletion move of a particle 'i' is accepted.
   This means removing the i'th row and column from the distance matrix.
---------------------------------------------------------------------------- */

void DistanceManager::update_remove(const unsigned int &i)
{
    remove_neigh(i);
    distance_mat.erase (distance_mat.begin() + i);
    distance_cube.erase (distance_cube.begin() + i);
    for (std::size_t j=0; j<distance_mat.size(); j++) {
        distance_mat[j].erase (distance_mat[j].begin() + i);
        distance_cube[j].erase (distance_cube[j].begin() + i);
    }
}


/* ----------------------------------------------------------------------------
   Update distance matrix after an insertion move of particle 'i'. This means
   appending a row and a column to the distance matrix. If more than one 
   particle is inserted, this function has to be called in ascending order
   for the particles.
---------------------------------------------------------------------------- */

void DistanceManager::update_insert(const unsigned int &i)
{
    double rij;
    std::valarray<double> delij;

    std::vector<double> rijs(box->npar);
    std::vector<std::valarray<double> > delijs(box->npar);

    // extend neighbor lists
    for (std::size_t k=0; k<ncutoff; k++) {
        //neigh_lists[k].insert(neigh_lists[k].begin() + i, {});
        neigh_lists[k].push_back({});
    }
    //for (std::vector<unsigned int> &neigh_list : neigh_lists) {
    //    neigh_list.push_back({});
    //}

    // extend distance matrix and update distance matrix and neighbor lists
    for (std::size_t j=0; j<box->npar-1; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        box->boundary->correct_distance(delij);
        rij = normsq(delij);
        delijs[j] = delij;
        rijs[j] = rij;
        distance_mat[j].insert(distance_mat[j].begin() + i, rij);
        distance_cube[j].insert(distance_cube[j].begin() + i, delij);
        if (i == j) continue;
        update_neigh(i, j, rij);
    }

    distance_mat.insert(distance_mat.begin() + i, rijs);
    distance_cube.insert(distance_cube.begin() + i, delijs);
}
