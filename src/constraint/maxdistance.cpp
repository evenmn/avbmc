/* ----------------------------------------------------------------------------
   Constraining the distance between two components A and B. All components A
   need to have a component B within the given radius in order to accept
   Monte Carlo moves. If the system is initialized such that this condition
   is not satisfied, no moves will be accepted (unless the constraint is 
   satisfied after a move
------------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>

#include "maxdistance.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../distance_manager.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Maximum distance constraint constructor. 'label1' and 'label2' are the IDs
   of the components A and B, respectively. 'rc_in' is the maximum distance
   between the components.
------------------------------------------------------------------------------- */

MaxDistance::MaxDistance(Box* box_in, std::string label1, std::string label2,
                         double rc_in) : Constraint(box_in)
{
    type1 = box->system->forcefield->label2type.at(label1);
    type2 = box->system->forcefield->label2type.at(label2);
    cutoff_id = box->distance_manager->add_cutoff(rc_in, label1, label2);
}


/* ----------------------------------------------------------------------------
   Verify that all particles of interest, 'particles',  meet the constraint
------------------------------------------------------------------------------- */

bool MaxDistance::verify()
{
    unsigned int i, typei;
    neigh_list = box->distance_manager->neigh_lists[cutoff_id];

    for (i=0; i<box->npar; i++) {
        typei = box->particles[i].type;
        if ((typei==type1 || typei==type2) && neigh_list[i].size()==0) {
            return false;
        }
    }
    return true;
}
