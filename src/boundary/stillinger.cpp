#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>
#include <memory>

#include "stillinger.h"
#include "../box.h"


/* ------------------------------------------------------
   Stillinger boundary constructor, initializing the
   Stillinger cluster criterion 'r_c'
--------------------------------------------------------- */

Stillinger::Stillinger(Box* box_in, double r_c_in)
    : Boundary(box_in)
{
    r_csq = r_c_in * r_c_in;
    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
    label = "Stillinger";
}


/*
Stillinger::Stillinger(std::shared_ptr<Box> box_in, double r_c_in)
    : Boundary(box_in)
{
    r_csq = r_c_in * r_c_in;
    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
    label = "Stillinger";
}
*/

/* ------------------------------------------------------
   Update neighbor lists of all particles according to
   the Stillinger criterion
--------------------------------------------------------- */

void Stillinger::update()
{
    neigh_lists.clear();
    for(int i=0; i<box->npar; i++){
        std::cout << "update" << i << std::endl;
        neigh_lists.push_back(box->build_neigh_list(i, r_csq));
    }
}


/* -------------------------------------------------------
   Check if particles are in cluster. 'checked' contains
   information about which atoms that we have checked 
   neighbor list of (to avoid circular check), and 
   'in_cluster' contains informations about which 
   particles that are part of cluster
---------------------------------------------------------- */

void Stillinger::check(const int i, std::valarray<int> &in_cluster, std::valarray<int> &checked)
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


/* ----------------------------------------------------------
   Check if Stillinger criterion is satisfied.
   If it is not satisfied, Monte Carlo moves
   should be rejected. Molecular dynamics
   simulations should abort in the same case.
------------------------------------------------------------- */

bool Stillinger::correct_position()
{
    std::cout << "correct_position1" << std::endl;
    std::cout << box->npar << std::endl;
    std::valarray<int> in_cluster2(0, 20);
    std::cout << "correct_position1" << std::endl;
    in_cluster.resize(box->npar, 0);
    std::cout << "correct_position1" << std::endl;
    checked.resize(box->npar, 0);
    std::cout << "correct_position2" << std::endl;
    //std::valarray<int> checked(0, box->npar);
    std::cout << "correct_position3" << std::endl;
    update();
    std::cout << "correct_position4" << std::endl;
    
    check(0, in_cluster, checked);
    std::cout << "correct_position5" << std::endl;

    if(in_cluster.sum() == box->npar){
        return true;
    }
    else{
        return false;
    }
    std::cout << "correct_position6" << std::endl;
}


/* ------------------------------------------------------------
   For Stillinger cluster criterion, the
   velocity does not need to be corrected
--------------------------------------------------------------- */

bool Stillinger::correct_velocity()
{
    return true;
}


/* ------------------------------------------------------------
   For Stillinger cluster criterion, the
   distance does not need to be corrected
--------------------------------------------------------------- */

bool Stillinger::correct_distance()
{
    return true;
}


/* ------------------------------------------------------------
   Get volume of cluster. This is done by first
   overestimating the volume as N * v_stillinger,
   and then subtracting the volume of overlapping
   sphere caps.
--------------------------------------------------------------- */
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
