#include "stillinger.h"
#include "../box.h"

Stillinger::Stillinger(class Box* box_in, double r_c_in)
    : Boundary(box_in)
{
    r_c = r_c_in;
    //v_c = 4 * datum::pi * pow(r_c, 3) / 3;
}

void Stillinger::update()
{
    /* Update neighbor lists of all particles
     */
    neigh_lists.clear();
    for(int i=0; i<box->npar; i++){
        neigh_lists.push_back(box->forcefield->build_neigh_list(i, r_c));
    }
}

void Stillinger::check(const int i, std::valarray<int> &in_cluster, std::valarray<int> &checked)
{
    /* Check if particles are in cluster. 'checked' contains
     * information about which atoms that we have checked 
     * neighbor list of (to avoid circular check), and 
     * 'in_cluster' contains informations about which 
     * particles that are part of cluster
     */
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


bool Stillinger::correct_position()
{
    /* Check if Stillinger criterion is satisfied.
     * If it is not satisfied, Monte Carlo moves
     * should be rejected. Molecular dynamics
     * simulations should abort in the same case.
     */

    std::valarray<int> in_cluster(0, box->npar);
    std::valarray<int> checked(0, box->npar);
    update();
    
    check(0, in_cluster, checked);

    if(in_cluster.sum() == box->npar){
        return true;
    }
    else{
        return false;
    }
}

bool Stillinger::correct_velocity()
{
    /* For Stillinger cluster criterion, the
     * velocity does not need to be corrected
     */
    return true;
}

bool Stillinger::correct_distance()
{
    /* For Stillinger cluster criterion, the
     * distance does not need to be corrected
     */
    return true;
}

/*
double Stillinger::comp_volume()
{
    // Get volume of cluster, but ignoring
    // the factor pi/3. This is done by first
    // overestimating the volume as N * v_stillinger,
    // and then subtracting the volume of overlapping
    // sphere caps.
    //

    double v_overest = box->npar * v_c;
    
    mat distance_mat_sqrd = box->forcefield->distance_mat;

    // compute height of all caps. This might be time-consuming
    mat h = (real(sqrtmat(distance_mat_sqrd)) / 2. - r_c);
    double v_caps = 2 * sum(sum(powmat(h, 2)%(-h + 3 * r_c)));

    return v_overest - v_caps;
}
*/
