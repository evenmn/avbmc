/* ----------------------------------------------------------------------------
   Constraining the number of neighbors component A has of component B.
------------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>

#include "minneigh.h"
#include "../box.h"
#include "../particle.h"
#include "../distance_manager.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Minimum neighbor constraint constructor. 'label1' and 'label2' are the IDs
   of the components A and B, respectively. 'rc_in' is the maximum number of
   neighbors that a particle of component A can have of particles of 
   component B.
------------------------------------------------------------------------------- */

MinNeigh::MinNeigh(Box* box_in, std::string label1, std::string label2,
                   double rc_in, int nc_in, bool restrict_lower_in)
    : Constraint(box_in)
{
    nc = nc_in;
    restrict_lower = restrict_lower_in;
    type1 = box->forcefield->label2type.at(label1);
    type2 = box->forcefield->label2type.at(label2);
    cutoff_id = box->distance_manager->add_cutoff(rc_in, label1, label2, true);
}


/* ----------------------------------------------------------------------------
   Verify that all particles of interest, 'particles',  meet the constraint
------------------------------------------------------------------------------- */

bool MinNeigh::verify()
{
    unsigned int i, typei, target_atoms;

    if (type1 == type2 && box->npartype[type1] < nc + 1) {
        // constraint does not apply if particle pair does not exist
        if (restrict_lower) {
            target_atoms = box->npartype[type1] - 1;
        }
        else {
            return true;
        }
    }
    else if (box->npartype[type2] < nc) {
        // type 1 cannot have nc neighbors of type2 if there is less than nc type2
        if (restrict_lower) {
            target_atoms = box->npartype[type2];
        }
        else {
            return true;
        }
    }
    else {
        target_atoms = nc;
    }
    /*
    if (type1 == type2 && box->npartype[type1] < nc + 1) {
        // constraint does not apply if particle pair does not exist
    }
    else if (box->npartype[type1] == 0 || box->npartype[type2] == 0) {
        // constraint does not apply if there is no particle of one type
    }
    else {
        neigh_list = box->distance_manager->neigh_lists[cutoff_id];
        for (i=0; i<box->npar; i++) {
            typei = box->particles[i].type;
            if (typei==type1  && neigh_list[i].size()<nc) {
                return false;
            }
        }
    }
    */
    neigh_list = box->distance_manager->neigh_lists[cutoff_id];
    for (i=0; i<box->npar; i++) {
        typei = box->particles[i].type;
        if (typei==type1  && neigh_list[i].size()<target_atoms) {
            return false;
        }
    }
    return true;
}
