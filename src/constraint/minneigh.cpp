/* ----------------------------------------------------------------------------
   Constraining the number of neighbors component A has of component B.
------------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>

#include "minneigh.h"
#include "../box.h"
#include "../system.h"
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
                   double rc_in, int nc_in) : Constraint(box_in)
{
    nc = nc_in;
    type1 = box->system->forcefield->label2type.at(label1);
    type2 = box->system->forcefield->label2type.at(label2);
    cutoff_id = box->distance_manager->add_cutoff(rc_in, label1, label2);
}


/* ----------------------------------------------------------------------------
   Verify that all particles of interest, 'particles',  meet the constraint
------------------------------------------------------------------------------- */

bool MinNeigh::verify()
{
    unsigned int i, typei;
    neigh_list = box->distance_manager->neigh_lists[cutoff_id];

    for (i=0; i<box->npar; i++) {
        typei = box->particles[i].type;
        if (typei==type1 && neigh_list[i].size()<nc) {
            return false;
        }
    }
    return true;
}
