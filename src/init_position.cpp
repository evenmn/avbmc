/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Initialize positions. Face-centered cube is currently the only crystal 
  initialization available. NB: A more powerful initialization method is to
  read an xyz-file, which can be implemented from an external software.
---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <valarray>
#include <string>

#include "init_position.h"
#include "particle.h"


/* ----------------------------------------------------------------------------
   Creating a face-centered cube of n^dim unit cells with 4 particles in each
   unit cell. The number of particles then becomes (dim+1) * n ^ dim. Each unit
   cell has a length d. L=nd
   Parameters
   ----------
   cells : int
       number of unit cells in each dimension
   lenbulk : float
       length of box
   dim : int
       number of dimensions
---------------------------------------------------------------------------- */

std::vector<std::valarray<double> > fcc(unsigned int ncells, double lenbulk,
    unsigned int ndim)
{
    unsigned int i, j, k;
    std::vector<std::valarray<double> > positions;

    if (ndim == 1) {
        for (i=0; i<ncells; i++) {
            positions.push_back({.0+i});
            positions.push_back({.5+i});
        }
    }
    else if (ndim == 2) {
        for (i=0; i<ncells; i++) {
            for (j=0; j<ncells; j++) {
                positions.push_back({.0+i, .0+j});
                positions.push_back({.0+i, .5+j});
                positions.push_back({.5+i, .0+j});
            }
        }
    }
    else if (ndim == 3) {
        for (i=0; i<ncells; i++) {
            for (j=0; j<ncells; j++) {
                for (k=0; k<ncells; k++) {
                    positions.push_back({.0+i, .0+j, .0+k});
                    positions.push_back({.0+i, .5+j, .5+k});
                    positions.push_back({.5+i, .0+j, .5+k});
                    positions.push_back({.5+i, .5+j, .0+k});
                }
            }
        }
    }
    else{
        std::cout << "The number of dimensions needs to be in (1,3)! "
                  << "Aborting" << std::endl;
        exit(0);
    }

    // scale properly
    double scale = lenbulk / ncells;
    for(std::valarray<double> &position : positions){
        position *= scale;
    }
    return positions;
}
