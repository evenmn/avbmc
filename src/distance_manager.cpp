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
------------------------------------------------------------------------------- */

DistanceManager::DistanceManager(Box* box_in, double cutoff_tol_in)
{
    ncutoff = 0;
    box = box_in;
    cutoff_tol = cutoff_tol_in;
}


/* ----------------------------------------------------------------------------
   Add a cutoff to the list of cutoffs in order to make the distance manager
   manage neighbor lists. This cutoff applies to all components. Returning the
   cutoff-ID, which has to be used when extracting the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double rc)
{
    double rcsq;
    unsigned int i, mode;

    mode = 0;
    rcsq = rc * rc;

    // check if a similar cutoff exists
    for (i=0; i<ncutoff; i++) {
        if (modes[i] == mode && fabs(cutoffs[i] - rcsq) < cutoff_tol) {
            return i;
        }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    modes.push_back(mode);
    cutoffs.push_back(rcsq);
    neigh_lists.push_back({});
    return ncutoff - 1;
}


/* ----------------------------------------------------------------------------
   Add a cutoff 'rc' to the list of cutoffs in order to make the distance
   manager manage neighbor lists. Here, the cutoff is restricted for two
   components 'label1' and 'label2'. If 'mutual' is true, both label1 and 
   label2's neighbor lists will be updated. Returning the cutoff-ID, which has
   to be used when extracting the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double rc, std::string label1,
                                         std::string label2, bool mutual)
{
    double rcsq;
    unsigned int i, j, type1, type2, mode;

    mode = 1;
    rcsq = rc * rc;
    type1 = box->forcefield->label2type.at(label1);
    type2 = box->forcefield->label2type.at(label2);

    // check if a similar cutoff exists
    j=0;
    for (i=0; i<ncutoff; i++) {
        if (modes[i] == mode) {
            if (types2[j] == type1 && types2[j] == type2 &&
                mutuals[j] == mutual && fabs(cutoffs[i] - rcsq) < cutoff_tol) {
                return i;
            }
            j++;
        }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    modes.push_back(mode);
    types1.push_back(type1);
    types2.push_back(type2);
    mutuals.push_back(mutual);
    cutoffs.push_back(rcsq);
    neigh_lists.push_back({});
    return ncutoff - 1;
}


/* ----------------------------------------------------------------------------
   Add a cutoff matrix 'rc' of shape (ntype, ntype) and create a combined
   neighbor list. Returning the cutoff-ID, which has to be used when extracting
   the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double **rc)
{
    unsigned int mode;

    mode = 2;
    ncutoff ++;
    modes.push_back(mode);
    cutoff_mats.push_back(rc);
    neigh_lists.push_back({});
    return ncutoff - 1;
}


/* ----------------------------------------------------------------------------
   Remove neighbor lists of and clear particle 'i' from the neighbor lists.
   This is done when a particle is moved or removed.
------------------------------------------------------------------------------- */

void DistanceManager::clear_neigh(unsigned int i)
{
    unsigned int k, l;
    std::vector<int>::iterator position;

    for (k=0; k<ncutoff; k++) {
        neigh_lists[k][i].clear();
        for (l=0; l<neigh_lists[k].size(); l++) {
            position = std::find(neigh_lists[k][l].begin(), neigh_lists[k][l].end(), i);
            if (position != neigh_lists[k][l].end()) {
                neigh_lists[k][l].erase(position);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Update neighbor lists of a particle pair 'i' and 'j'. This is done when a
   particle is moved or inserted.
------------------------------------------------------------------------------- */

void DistanceManager::update_neigh(unsigned int i, unsigned int j, double rij)
{
    unsigned int k, l, m, typei, typej;

    typei = box->particles[i].type;
    typej = box->particles[j].type;

    l=m=0;
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
}


/* ----------------------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
------------------------------------------------------------------------------- */

double DistanceManager::normsq(std::valarray<double> array)
{
    double norm;
    unsigned int i;

    norm = 0.;
    for (i=0; i < array.size(); i++)
    {
        norm += array[i] * array[i];
    }
    return norm;
}


/* ----------------------------------------------------------------------------
   Initialize distance matrix and neighbor lists
------------------------------------------------------------------------------- */

void DistanceManager::initialize()
{
    double rij;
    unsigned int i, j, k, npar;
    std::valarray<double> posi, delij;

    // initialize matrices
    npar = box->npar;
    distance_mat.resize(npar);
    distance_cube.resize(npar);
    for (i=0; i<npar; i++) {
        distance_mat[i].resize(npar);
        distance_cube[i].resize(npar);
        for (k=0; k<ncutoff; k++) {
            neigh_lists[k].push_back({});
        }
    }
    
    // fill matrices
    for (i=0; i<npar; i++) {
        posi = box->particles[i].r;
        for (j=0; j<i; j++) {
            delij = box->particles[j].r - box->particles[i].r;
            box->boundary->correct_distance(delij);
            rij = normsq(delij);
            distance_mat[i][j] = rij;
            distance_cube[i][j] = delij;
            distance_mat[j][i] = rij;
            distance_cube[j][i] = -delij;
            update_neigh(i, j, rij);
            //update_neigh(j, i, rij);
        }
    }
}


/* ----------------------------------------------------------------------------
   Store old matrices
------------------------------------------------------------------------------- */

void DistanceManager::set()
{
    distance_mat_old = distance_mat;
    distance_cube_old = distance_cube;
    neigh_lists_old = neigh_lists;
}


/* ----------------------------------------------------------------------------
   Reset new matrices
------------------------------------------------------------------------------- */

void DistanceManager::reset()
{
    distance_mat = distance_mat_old;
    distance_cube = distance_cube_old;
    neigh_lists = neigh_lists_old;
}


/* ----------------------------------------------------------------------------
   Update distance matrix due to translational move of a particle 'i'. This 
   corresponds to updating the i'th row and column in the distance matrix.
------------------------------------------------------------------------------- */

void DistanceManager::update_trans(unsigned int i)
{
    double rij;
    unsigned int j;
    std::valarray<double> delij;

    clear_neigh(i);
    for (j=0; j<i; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        box->boundary->correct_distance(delij);
        rij = normsq(delij);
        distance_mat[i][j] = rij;
        distance_cube[i][j] = delij;
        distance_mat[j][i] = rij;
        distance_cube[j][i] = -delij;
        update_neigh(i, j, rij);
    }
    for (j=i+1; j<box->npar; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        box->boundary->correct_distance(delij);
        rij = normsq(delij);
        distance_mat[i][j] = rij;
        distance_cube[i][j] = delij;
        distance_mat[j][i] = rij;
        distance_cube[j][i] = -delij;
        update_neigh(i, j, rij);
    }
}


/* ----------------------------------------------------------------------------
   Update distance matrix after a deletion move of a particle 'i' is accepted.
   This means removing the i'th row and column from the distance matrix.
------------------------------------------------------------------------------- */

void DistanceManager::update_remove(unsigned int i)
{
    unsigned int j;

    clear_neigh(i);
    distance_mat.erase (distance_mat.begin() + i);
    distance_cube.erase (distance_cube.begin() + i);
    for (j=0; j<distance_mat.size(); j++) {
        distance_mat[j].erase (distance_mat[j].begin() + i);
        distance_cube[j].erase (distance_cube[j].begin() + i);
    }
    std::cout << "NEW SIZE: " << distance_mat.size() << "x" << distance_mat[0].size() << std::endl;
}


/* ----------------------------------------------------------------------------
   Update distance matrix after an insertion move of particle 'i'. This means
   appending a row and a column to the distance matrix. If more than one 
   particle is inserted, this function has to be called in ascending order
   for the particles.
------------------------------------------------------------------------------- */

void DistanceManager::update_insert(unsigned int i)
{
    double rij;
    unsigned int j, k;
    std::valarray<double> delij;

    std::vector<double> rijs(box->npar);
    std::vector<std::valarray<double> > delijs(box->npar);

    // extend neighbor lists
    for (k=0; k<ncutoff; k++) {
        //neigh_lists[k].insert(neigh_lists[k].begin() + i, {});
        neigh_lists[k].push_back({});
    }

    // extend distance matrix and update distance matrix and neighbor lists
    for (j=0; j<box->npar-1; j++) {
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
