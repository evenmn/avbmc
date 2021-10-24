#include "init_position.h"


mat fcc(int ncells, double lenbulk, int ndim)
{
    /* Creating a face-centered cube of n^dim unit cells with
     * 4 particles in each unit cell. The number of particles
     * then becomes (dim+1) * n ^ dim. Each unit cell has a
     * length d. L=nd
     * Parameters
     * ----------
     * cells : int
     *     number of unit cells in each dimension
     * lenbulk : float
     *     length of box
     * dim : int
     *     number of dimensions
     */
    int npar = (ndim+1) * pow(ncells, ndim);
    mat r = zeros(npar, ndim);
    int counter = 0;
    if(ndim == 1){
        for(int i=0; i<ncells; i++){
            r.row(counter+0) = {.0+i};
            r.row(counter+1) = {.5+i};
            counter += 2;
        }
    }
    else if(ndim == 2){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                r.row(counter+0) = {.0+i, .0+j};
                r.row(counter+1) = {.0+i, .5+j};
                r.row(counter+2) = {.5+i, .0+j};
                counter += 4;
            }
        }
    }
    else if(ndim == 3){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                for(int k=0; k<ncells; k++){
                    r.row(counter+0) = {.0+i, .0+j, .0+k};
                    r.row(counter+1) = {.0+i, .5+j, .5+k};
                    r.row(counter+2) = {.5+i, .0+j, .5+k};
                    r.row(counter+3) = {.5+i, .5+j, .0+k};
                    counter += 4;
                }
            }
        }
    }
    else{
        cout << "The number of dimensions needs to be in (1,3), aborting" << endl;
        exit(0);
    }

    // scale initial positions correctly
    r *= lenbulk / ncells;
    return r;
}

mat from_xyz(string filename)
{
    mat r(1, 3);
    return r;
}
