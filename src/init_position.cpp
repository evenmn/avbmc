#include <iostream>
#include <vector>
#include <valarray>
#include <string>

#include "init_position.h"
#include "particle.h"


/* -----------------------------------------------------------
   Creating a face-centered cube of n^dim unit cells with
   4 particles in each unit cell. The number of particles
   then becomes (dim+1) * n ^ dim. Each unit cell has a
   length d. L=nd
   Parameters
   ----------
   cells : int
       number of unit cells in each dimension
   lenbulk : float
       length of box
   dim : int
       number of dimensions
--------------------------------------------------------------- */

std::vector<std::valarray<double> > fcc(int ncells, double lenbulk, int ndim)
{
    // initialize position vector
    std::vector<std::valarray<double> > positions;

    if(ndim == 1){
        for(int i=0; i<ncells; i++){
            positions.push_back({.0+i});
            positions.push_back({.5+i});
        }
    }
    else if(ndim == 2){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                positions.push_back({.0+i, .0+j});
                positions.push_back({.0+i, .5+j});
                positions.push_back({.5+i, .0+j});
            }
        }
    }
    else if(ndim == 3){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                for(int k=0; k<ncells; k++){
                    positions.push_back({.0+i, .0+j, .0+k});
                    positions.push_back({.0+i, .5+j, .5+k});
                    positions.push_back({.5+i, .0+j, .5+k});
                    positions.push_back({.5+i, .5+j, .0+k});
                }
            }
        }
    }
    else{
        std::cout << "The number of dimensions needs to be in (1,3), aborting" << std::endl;
        exit(0);
    }

    // scale properly
    double scale = lenbulk / ncells;
    for(std::valarray<double> &position : positions){
        position *= scale;
    }

    return positions;
}

std::vector<class Particle *> from_xyz(std::string filename)
{
    std::vector<class Particle *> particles;  // = read_xyz(filename);
    return particles;
}
