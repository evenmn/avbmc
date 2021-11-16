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
    for(int i=0; i<box->npar; i++){
        cout << "i: " << i << endl;
        neigh_lists.push_back(box->forcefield->build_neigh_list(i, r_c));
    }
}

void Stillinger::check(const int i, std::valarray<bool> &in_cluster, std::valarray<bool> &checked)
{
    /* Check if particles are in cluster. 'checked' contains
     * information about which atoms that we have checked 
     * neighbor list of (to avoid circular check), and 
     * 'in_cluster' contains informations about which 
     * particles that are part of cluster
     */
    cout << i << endl;
    cout << in_cluster.size() << endl;
    cout << checked.size() << endl;
    if(!checked[i]){
        checked[i] = true;
        in_cluster[i] = true;
        for(int j : neigh_lists[i]){
            cout << "j: " << j << endl;
            in_cluster[j] = true;
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

    std::valarray<bool> in_cluster(false, box->npar);
    std::valarray<bool> checked(false, box->npar);
    cout << "correct_position1" << endl;
    update();
    
    cout << "correct_position2" << endl;
    check(0, in_cluster, checked);
    cout << "correct_position3" << endl;
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
