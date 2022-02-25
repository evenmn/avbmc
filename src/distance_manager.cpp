#include <iostream>
#include <vector>
#include <valarray>
#include <algorithm>  // for removing vector element by value

#include "distance_manager.h"
#include "system.h"
#include "box.h"


/* ----------------------------------------------------------------------------
   Distance manager class constructor. 
------------------------------------------------------------------------------- */

DistanceManager::DistanceManager(Box* box_in, double cutoff_tol_in)
{
    ncutoff = 0;
    box = box_in;
    cutoff_tol = cutoff_tol_in;
}


/* ----------------------------------------------------------------------------
   Add a cutoff to the list of cutoffs in order for the distance manager
   to make neighbor lists. Returning the cutoff-ID, which has to be used
   when extracting the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double rc)
{
    double rcsq;
    bool existing_cutoff;
    unsigned int i, cutoffid;

    rcsq = rc * rc;

    // check if a similar cutoff exists
    existing_cutoff = false;
    for (i=0; i<ncutoff; i++) {
        if (fabs(cutoffs[i] - rcsq) < cutoff_tol) {
            existing_cutoff = true;
            cutoffid = i;
            break;
        }
    }

    // add cutoff to list of cutoffs if it does
    // not already exist
    if (!existing_cutoff) {
        cutoffid = ncutoff;
        ncutoff ++;
        cutoffs.push_back(rcsq);
        neigh_lists.push_back({});
    }
    return cutoffid;
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
   Compute distance between two particles 'i' and 'j' and update distance
   matrix and neighbor lists for these two specific particles
------------------------------------------------------------------------------- */
/*
void DistanceManager::update(int i, int j)
{
    double rij, cutoff;
    valarray<double> delij;

    delij = box->particles[i].r - box->particles[j].r;
    rij = normsq(delij);
    distance_mat[i][j] = rij;
    distance_cube[i][j] = delij;
    for (cutoff : cutoffs) {
        if (rij < cutoff) {
            neigh_lists[i].push_back(j);
        }
    }
}
*/

/* ----------------------------------------------------------------------------
   Initialize distance matrix and neighbor lists
------------------------------------------------------------------------------- */

void DistanceManager::initialize()
{
    std::cout << "initialize1" << std::endl;
    double rij;
    unsigned int i, j, k, npar;
    std::valarray<double> posi, delij;

    // initialize matrices
    npar = box->npar;
    distance_mat.resize(npar);
    distance_cube.resize(npar);
    for (k=0; k<npar; k++) {
        distance_mat[k].resize(npar);
        distance_cube[k].resize(npar);
    }
    
    // fill matrices
    for (i=0; i<npar; i++) {
        posi = box->particles[i].r;
        for (j=0; j<i; j++) {
            delij = box->particles[j].r - box->particles[i].r;
            rij = normsq(delij);
            distance_mat[i][j] = rij;
            distance_cube[i][j] = delij;
            distance_mat[j][i] = rij;
            distance_cube[j][i] = -delij;
            for (k=0; k<ncutoff; k++) {
                if (rij < cutoffs[k]) {
                    neigh_lists[k][i].push_back(j);
                    neigh_lists[k][j].push_back(i);
                }
            }
        }
    }
    std::cout << "initialize2" << std::endl;
}


/* ----------------------------------------------------------------------------
   Update distance matrix due to translational move of a particle 'i'. This 
   corresponds to updating the i'th row and column in the distance matrix.
------------------------------------------------------------------------------- */

void DistanceManager::update_trans(unsigned int i)
{
    double rij;
    unsigned int j, k, l;
    std::valarray<double> delij;

    // clear neighbor lists of particle i and
    // remove particle i from all neigh lists
    for (k=0; k<ncutoff; k++) {
        neigh_lists[k][i].clear();
        for (l=0; l<neigh_lists[k].size(); l++) {
            std::vector<int> nl = neigh_lists[k][l];
            std::vector<int>::iterator position = std::find(nl.begin(), nl.end(), i);
            if (position != nl.end()) {
                neigh_lists[k][l].erase(position);
            }
        }
    }

    // update distance matrix and neighbor lists
    for (j=0; j<box->npar; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        rij = normsq(delij);
        distance_mat[i][j] = rij;
        distance_cube[i][j] = delij;
        distance_mat[j][i] = rij;
        distance_cube[j][i] = -delij;
        for (k=0; k<ncutoff; k++) {
            if (rij < cutoffs[k]) {
                neigh_lists[k][i].push_back(j);
                neigh_lists[k][j].push_back(i);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Update distance matrix after a deletion move of a particle 'i' is accepted.
   This means removing the i'th row and column from the distance matrix.
------------------------------------------------------------------------------- */

void DistanceManager::update_remove(unsigned int i)
{
    unsigned int j, k, l;

    // remove neighbor lists of particle i and remove
    // i from all neighbor lists
    for (k=0; k<ncutoff; k++) {
        neigh_lists[k].erase(neigh_lists[k].begin() + i);
        for (l=0; l<neigh_lists[k].size(); l++) {
            std::vector<int> nl = neigh_lists[k][l];
            std::vector<int>::iterator position = std::find(nl.begin(), nl.end(), i);
            if (position != nl.end()) {
                neigh_lists[k][l].erase(position);
            }
        }
    }

    // remove row and column from distance matrix
    distance_mat.erase (distance_mat.begin() + i);
    distance_cube.erase (distance_cube.begin() + i);
    for (j=0; j<distance_mat.size(); j++) {
        distance_mat[j].erase (distance_mat[j].begin() + i);
        distance_cube[j].erase (distance_cube[j].begin() + i);
    }
}


/* ----------------------------------------------------------------------------
   Update distance matrix after an insertion move of particle 'i'. This means
   appending a row and a column to the distance matrix. If more than one 
   particle is inserted, this function has to be called in ascending order
   for the particles.
------------------------------------------------------------------------------- */

void DistanceManager::update_insert(unsigned int i)
{
    double rij, cutoff;
    unsigned int j, k;
    std::valarray<double> delij;

    std::vector<double> rijs(box->npar);
    std::vector<std::valarray<double> > delijs(box->npar);

    // extend neighbor lists
    for (k=0; k<ncutoff; k++) {
        neigh_lists[k].insert(neigh_lists[k].begin() + i, {});
    }

    // extend distance matrix and update distance matrix and neighbor lists
    for (j=0; j<box->npar; j++) {
        delij = box->particles[j].r - box->particles[i].r;
        rij = normsq(delij);
        delijs[j] = delij;
        rijs[j] = rij;
        distance_mat[j].insert(distance_mat[j].begin() + i, rij);
        distance_cube[j].insert(distance_cube[j].begin() + i, delij);
        for (k=0; k<ncutoff; k++) {
            if (rij < cutoffs[k]) {
                neigh_lists[k][i].push_back(j);
                neigh_lists[k][j].push_back(i);
            }
        }
    }
    distance_mat.insert(distance_mat.begin() + i, rijs);
    distance_cube.insert(distance_cube.begin() + i, delijs);
}
