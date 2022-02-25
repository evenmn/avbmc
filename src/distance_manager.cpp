#include <iostream>
#include <vector>
#include <valarray>
#include <algorithm>  // for removing vector element by value

#include "distance_manager.h"
#include "box.h"
#include "system.h"
#include "forcefield/forcefield.h"


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
   Add a cutoff to the list of cutoffs in order to make the distance manager
   manage neighbor lists. This cutoff applies to all components. Returning the
   cutoff-ID, which has to be used when extracting the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double rc)
{
    double rcsq;
    unsigned int i;

    rcsq = rc * rc;

    // check if a similar cutoff exists
    for (i=0; i<ncutoff; i++) {
        if (!res[i] && fabs(cutoffs[i] - rcsq) < cutoff_tol) {
            return i;
        }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    res.push_back(false);
    cutoffs.push_back(rcsq);
    neigh_lists.push_back({});
    return ncutoff - 1;
}


/* ----------------------------------------------------------------------------
   Add a cutoff 'rc' to the list of cutoffs in order to make the distance
   manager manage neighbor lists. Here, the cutoff is restricted for two
   components 'label1' and 'label2'. Returning the cutoff-ID, which has to be
   used when extracting the correct neighbor list.
------------------------------------------------------------------------------- */

unsigned int DistanceManager::add_cutoff(double rc, std::string label1,
                                         std::string label2)
{
    double rcsq;
    unsigned int i, j, type1, type2;

    rcsq = rc * rc;
    type1 = box->system->forcefield->label2type.at(label1);
    type2 = box->system->forcefield->label2type.at(label2);

    // check if a similar cutoff exists
    j=0;
    for (i=0; i<ncutoff; i++) {
        if (res[i]) {
            if (types2[j] == type1 && types2[j] == type2 &&
                fabs(cutoffs[i] - rcsq) < cutoff_tol) {
                return i;
            }
            j ++;
        }
    }

    // add cutoff to list of cutoffs if it does not already exist
    ncutoff ++;
    res.push_back(true);
    types1.push_back(type1);
    types2.push_back(type2);
    cutoffs.push_back(rcsq);
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
    std::vector<int> nl;
    std::vector<int>::iterator position;

    for (k=0; k<ncutoff; k++) {
        neigh_lists[k][i].clear();
        for (l=0; l<neigh_lists[k].size(); l++) {
            nl = neigh_lists[k][l];
            position = std::find(nl.begin(), nl.end(), i);
            if (position != nl.end()) {
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
    unsigned int k, l, typei, typej;

    typei = box->particles[i].type;
    typej = box->particles[j].type;

    l=0;
    for (k=0; k<ncutoff; k++) {
        if (res[k]) { 
            if (types1[l]==typei && types2[l]==typej
                && rij < cutoffs[k]) {
                neigh_lists[k][i].push_back(j);
                neigh_lists[k][j].push_back(i);
            }
            l++;
        }
        else {
            if (rij < cutoffs[k]) {
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
    unsigned int i, j, k, l, typei, typej, npar;
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
            typej = box->particles[j].type;
            delij = box->particles[j].r - box->particles[i].r;
            rij = normsq(delij);
            distance_mat[i][j] = rij;
            distance_cube[i][j] = delij;
            distance_mat[j][i] = rij;
            distance_cube[j][i] = -delij;
            update_neigh(i, j, rij);
        }
    }
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
    for (j=0; j<box->npar; j++) {
        delij = box->particles[j].r - box->particles[i].r;
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
        update_neigh(i, j, rij);
    }
    distance_mat.insert(distance_mat.begin() + i, rijs);
    distance_cube.insert(distance_cube.begin() + i, delijs);
}
