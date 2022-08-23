/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-23)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Constraining the number of neighbors component A has of component B.
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>

#include "maxneigh.h"
#include "../box.h"
#include "../particle.h"
#include "../distance_manager.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Maximum neighbor constraint constructor. 'label1' and 'label2' are the IDs
   of the components A and B, respectively. 'rc_in' is the maximum number of
   neighbors that a particle of component A can have of particles of 
   component B.
---------------------------------------------------------------------------- */

MaxNeigh::MaxNeigh(Box* box_in, std::string label1, std::string label2,
                   double rc_in, int nc_in) : Constraint(box_in)
{
    nc = nc_in;
    type1 = box->forcefield->label2type.at(label1);
    type2 = box->forcefield->label2type.at(label2);
    cutoff_id = box->distance_manager->add_cutoff(rc_in, label1, label2);
}


/* ----------------------------------------------------------------------------
   Verify that all particles of interest, 'particles',  meet the constraint
---------------------------------------------------------------------------- */

bool MaxNeigh::verify()
{
    unsigned int i, typei;
    auto t0 = Time::now();

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
            if (typei==type1 && neigh_list[i].size()>nc) {
                nreject ++;
                auto t1 = Time::now();
                fsec fs = t1 - t0;
                cum_time += fs.count();
                return false;
            }
        }
    }
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    cum_time += fs.count();
    return true;
}
