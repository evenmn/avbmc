#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "stillinger.h"
#include "../forcefield/forcefield.h"
#include "../box.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   Stillinger boundary constructor, initializing the Stillinger cluster
   criterion 'r_c'
------------------------------------------------------------------------------- */

Stillinger::Stillinger(Box* box_in, double r_c_in)
    : Boundary(box_in)
{
    r_csq = r_c_in * r_c_in;
    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
    label = "Stillinger of radius " + std::to_string(r_c_in);

    // fill r_csq_mat with r_csq
    ntype = box->system->forcefield->ntype;
    r_csq_mat = new double*[ntype];
    for (unsigned int i=0; i<ntype; i++) {
        r_csq_mat[i] = new double[ntype];
        for (unsigned int j=0; j<ntype; j++) {
            r_csq_mat[i][j] = r_csq;
        }
    }
}


/* ----------------------------------------------------------------------------
   Set stillinger criterion for two substances 'label1' and 'label2'. This
   will overwrite the critical distance given when constructing the object.
------------------------------------------------------------------------------- */

void Stillinger::set_crit(std::string label1, std::string label2, double r_c)
{
    unsigned int type1 = box->system->forcefield->label2type[label1];
    unsigned int type2 = box->system->forcefield->label2type[label2];
    r_csq_mat[type1][type2] = r_c * r_c;
    r_csq_mat[type2][type1] = r_c * r_c;
}


/* ----------------------------------------------------------------------------
   Update neighbor lists of all particles according to the Stillinger criterion
------------------------------------------------------------------------------- */

void Stillinger::update()
{
    neigh_lists.clear();
    for(unsigned int i=0; i<box->npar; i++){
        neigh_lists.push_back(box->build_neigh_list(i, r_csq_mat));
    }
}


/* ----------------------------------------------------------------------------
   Check if particles are in cluster. 'checked' contains information about
   which atoms that we have checked neighbor list of (to avoid circular check),
   and 'in_cluster' contains informations about which particles that are part
   of cluster
------------------------------------------------------------------------------- */

void Stillinger::check(const int i, std::valarray<int> &in_cluster,
                       std::valarray<int> &checked)
{
    if(!checked[i]){
        checked[i] = 1;
        in_cluster[i] = 1;
        for(int j : neigh_lists[i]){
            in_cluster[j] = 1;
        }
        for(int j : neigh_lists[i]){
            check(j, in_cluster, checked);
        }
    }
}


/* ----------------------------------------------------------------------------
   Check if Stillinger criterion is satisfied. If it is notsatisfied, Monte
   Carlo moves should be rejected. Molecular dynamics simulations should abort
   in the same case.
------------------------------------------------------------------------------- */

bool Stillinger::correct_position()
{
    in_cluster.resize(box->npar, 0);
    checked.resize(box->npar, 0);
    update();
    
    check(0, in_cluster, checked);

    if(in_cluster.sum() == box->npar){
        return true;
    }
    else{
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

Stillinger::~Stillinger()
{
    for (unsigned int i=0; i<ntype; i++) {
        delete[] r_csq_mat[i];
    }
    delete[] r_csq_mat;
}
