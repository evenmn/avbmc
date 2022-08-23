/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-23)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Cluster criterion introduced by F. H. Stillinger in his "Rigorous Basis of
   the Frenkel‐Band Theory of Association Equilibrium" (1963). This constraint
   forces all particles to be part of the same cluster by rejecting all moves
   that splits up the cluster. Note that if the system is initialized with 
   more than one cluster, all moves will be rejected unless a particle is 
   inserted such that all clusters are bridged. 
---------------------------------------------------------------------------- */

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
---------------------------------------------------------------------------- */

Stillinger::Stillinger(Box* box_in, double rc)
    : Constraint(box_in)
{
    unsigned int i, j;

    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
    v_c = 0.;
    ntype = 0;
    label = "Stillinger";

    // fill r_csq_mat with r_csq
    ntype = box->forcefield->ntype;
    r_csq_mat = new double*[ntype];
    for (i=0; i<ntype; i++) {
        r_csq_mat[i] = new double[ntype];
        for (j=0; j<ntype; j++) {
            r_csq_mat[i][j] = rc * rc;
        }
    }
    cutoff_id = box->distance_manager->add_cutoff(r_csq_mat);
    vecid = box->distance_manager->mapid2vector[cutoff_id];
}


/* ----------------------------------------------------------------------------
   Copy constructor, needed to fullfil the rule of three, Marshall Cline (1991)
---------------------------------------------------------------------------- */

Stillinger::Stillinger(const Stillinger &other) 
    : Constraint(other), neigh_lists(other.neigh_lists),
      in_cluster(other.in_cluster), checked(other.checked)
{
    unsigned int i, j;

    ntype = other.ntype;
    v_c = other.v_c;
    r_csq_mat = new double*[ntype];
    for (i=0; i<ntype; i++) {
        r_csq_mat[i] = new double[ntype];
        for (j=0; j<ntype; j++) {
            r_csq_mat[i][j] = other.r_csq_mat[i][j];
        }
    } 
}


/* ----------------------------------------------------------------------------
   Swap function, used by the copy assignment operator
---------------------------------------------------------------------------- */

void Stillinger::swap(Stillinger &other)
{
    unsigned int ntype_tmp = ntype;
    ntype = other.ntype;
    other.ntype = ntype_tmp;
    unsigned int cutoff_id_tmp = cutoff_id;
    cutoff_id = other.cutoff_id;
    other.cutoff_id = cutoff_id_tmp;
    double v_c_tmp = v_c;
    v_c = other.v_c;
    other.v_c = v_c_tmp;
    double **r_csq_mat_tmp = r_csq_mat;
    r_csq_mat = other.r_csq_mat;
    other.r_csq_mat = r_csq_mat_tmp;
    std::vector<std::vector<int> > neigh_lists_tmp = neigh_lists;
    neigh_lists = other.neigh_lists;
    other.neigh_lists = neigh_lists_tmp;
    std::valarray<char> in_cluster_tmp = in_cluster;
    in_cluster = other.in_cluster;
    other.in_cluster = in_cluster_tmp;
    std::valarray<char> checked_tmp = checked;
    checked = other.checked;
    other.checked = checked_tmp;
}


/* ----------------------------------------------------------------------------
   Set stillinger criterion for two substances 'label1' and 'label2'. This
   will overwrite the critical distance given when constructing the object.
---------------------------------------------------------------------------- */

void Stillinger::set_criterion(std::string label1, std::string label2, double rc)
{
    unsigned int type1 = box->forcefield->label2type[label1];
    unsigned int type2 = box->forcefield->label2type[label2];
    
    box->distance_manager->cutoff_mats[vecid][type1][type2] = rc * rc;
    box->distance_manager->cutoff_mats[vecid][type2][type1] = rc * rc;
}


/* ----------------------------------------------------------------------------
   Check if particles are in cluster. 'checked' contains information about
   which atoms that we have checked neighbor list of (to avoid circular check),
   and 'in_cluster' contains informations about which particles that are part
   of cluster.
---------------------------------------------------------------------------- */

void Stillinger::check_neigh_recu(const int i, std::valarray<char> &in_cluster,
                                  std::valarray<char> &checked)
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
---------------------------------------------------------------------------- */

bool Stillinger::verify()
{
    auto t0 = Time::now();
    in_cluster.resize(box->npar, 0);
    checked.resize(box->npar, 0);
    
    check_neigh_recu(0, in_cluster, checked);
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    cum_time += fs.count();

    if ((unsigned int) in_cluster.sum() == box->npar) {
        return true;
    }
    else {
        nreject ++;
        return false;
    }
}


/* ----------------------------------------------------------------------------
   Get volume of cluster. This is done by first overestimating the volume as
   N * v_stillinger, and then subtracting the volume of overlapping sphere
   caps.
---------------------------------------------------------------------------- */
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
---------------------------------------------------------------------------- */

Stillinger::~Stillinger()
{
    unsigned int i;

    for (i=0; i<ntype; i++) {
        delete[] r_csq_mat[i];
    }
    delete[] r_csq_mat;
}
