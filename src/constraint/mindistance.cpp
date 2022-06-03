/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

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

#include "mindistance.h"
#include "../box.h"
#include "../particle.h"
#include "../distance_manager.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Minimum distance constraint constructor. 'label1' and 'label2' are the IDs
   of the components A and B, respectively. 'rc_in' is the minimum distance
   between the components.
------------------------------------------------------------------------------- */

MinDistance::MinDistance(Box* box_in, std::string label1, std::string label2,
                         double rc_in) : Constraint(box_in)
{
    type1 = box->forcefield->label2type.at(label1);
    type2 = box->forcefield->label2type.at(label2);
    cutoff_id = box->distance_manager->add_cutoff(rc_in, label1, label2);
}


/* ----------------------------------------------------------------------------
   Verify that all particles of interest, 'particles',  meet the constraint
------------------------------------------------------------------------------- */

bool MinDistance::verify()
{
    unsigned int i, typei;

    if (type1 == type2 && box->npartype[type1] < 2) {
        // constraint does not apply if particle pair does not exist
    }
    else if (box->npartype[type1] == 0 || box->npartype[type2] == 0) {
        // constraint does not apply if there is no particle of one type
    }
    else {
        neigh_list = box->distance_manager->neigh_lists[cutoff_id];

        for (i=0; i<box->npar; i++) {
            typei = box->particles[i].type;
            if ((typei==type1 || typei==type2) && neigh_list[i].size()>0) {
                return false;
            }
        }
    }
    return true;
}
