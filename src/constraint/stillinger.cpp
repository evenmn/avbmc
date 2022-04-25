/* ----------------------------------------------------------------------------
   Cluster criterion introduced by F. H. Stillinger in his "Rigorous Basis of
   the Frenkel‐Band Theory of Association Equilibrium" (1963). This constraint
   forces all particles to be part of the same cluster by rejecting all moves
   that splits up the cluster. Note that if the system is initialized with 
   more than one cluster, all moves will be rejected unless a particle is 
   inserted such that all clusters are bridged. 
------------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "stillinger.h"
#include "../box.h"
#include "../distance_manager.h"
#include "../forcefield/forcefield.h"
#include "../constraint/constraint.h"


/* ----------------------------------------------------------------------------
   Stillinger criterion constructor, initializing the Stillinger cluster
   criterion 'rc'. This contraint differs from the 'max distance' constraint
   in the way that it ensures that the system consists of exactly one cluster
   at all times. 
------------------------------------------------------------------------------- */

Stillinger::Stillinger(Box* box_in, double rc)
    : Constraint(box_in)
{
    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
    label = "Stillinger of radius " + std::to_string(rc);

    // fill r_csq_mat with r_csq
    ntype = box->forcefield->ntype;
    r_csq_mat = new double*[ntype];
    for (unsigned int i=0; i<ntype; i++) {
        r_csq_mat[i] = new double[ntype];
        for (unsigned int j=0; j<ntype; j++) {
            r_csq_mat[i][j] = rc * rc;
        }
    }
    cutoff_id = box->distance_manager->add_cutoff(r_csq_mat);
}


/* ----------------------------------------------------------------------------
   Set stillinger criterion for two substances 'label1' and 'label2'. This
   will overwrite the critical distance given when constructing the object.
------------------------------------------------------------------------------- */

void Stillinger::set_criterion(std::string label1, std::string label2, double rc)
{
    unsigned int type1 = box->forcefield->label2type[label1];
    unsigned int type2 = box->forcefield->label2type[label2];
    
    box->distance_manager->cutoff_mats[cutoff_id][type1][type2] = rc * rc;
    box->distance_manager->cutoff_mats[cutoff_id][type2][type1] = rc * rc;
}


/* ----------------------------------------------------------------------------
   Check if particles are in cluster. 'checked' contains information about
   which atoms that we have checked neighbor list of (to avoid circular check),
   and 'in_cluster' contains informations about which particles that are part
   of cluster.
------------------------------------------------------------------------------- */

void Stillinger::check_neigh_recu(const int i, std::valarray<int> &in_cluster,
                                  std::valarray<int> &checked)
{
    if (!checked[i]) {
        checked[i] = 1;
        in_cluster[i] = 1;
        for (int j : box->distance_manager->neigh_lists[cutoff_id][i]) {
            in_cluster[j] = 1;
        }
        for (int j : box->distance_manager->neigh_lists[cutoff_id][i]) {
            check_neigh_recu(j, in_cluster, checked);
        }
    }
}


/* ----------------------------------------------------------------------------
   Check if Stillinger criterion is satisfied. If it is not satisfied, Monte
   Carlo moves should be rejected. Molecular dynamics simulations should abort
   in the same case.
------------------------------------------------------------------------------- */

bool Stillinger::verify()
{
    in_cluster.resize(box->npar, 0);
    checked.resize(box->npar, 0);
    
    check_neigh_recu(0, in_cluster, checked);

    if (in_cluster.sum() == box->npar) {
        return true;
    }
    else {
        return false;
    }
}


/* ----------------------------------------------------------------------------
   Get volume of cluster. This is done by first overestimating the volume as
   N * v_stillinger, and then subtracting the volume of overlapping sphere
   caps.
------------------------------------------------------------------------------- */
/*
double Stillinger::comp_volume()
{
    double v_overest = box->npar * v_c;
    
    mat distance_mat_sqrd = box->forcefield->distance_mat;

    // compute height of all caps. This might be time-consuming
    mat h = (real(sqrtmat(distance_mat_sqrd)) / 2. - r_c);
    double v_caps = 2 * sum(sum(powmat(h, 2)%(-h + 3 * r_c)));

    return v_overest - v_caps;
}
*/

/* ----------------------------------------------------------------------------
   Stillinger destructor, frees up allocated memory
------------------------------------------------------------------------------- */

Stillinger::~Stillinger()
{
    unsigned int i;
    for (i=0; i<ntype; i++) {
        delete[] r_csq_mat[i];
    }
    delete[] r_csq_mat;
}
