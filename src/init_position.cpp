#include "init_position.h"
#include "particle.h"


std::vector<Particle *> fcc(int ncells, double lenbulk, int ndim)
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

    // initialize particle vector
    int npar = (ndim+1) * pow(ncells, ndim);
    std::vector<Particle *> particles;
    for(int i=0; i<npar; i++){
        Particle *particle;
        particles.push_back(particle);
    }

    int counter = 0;
    if(ndim == 1){
        for(int i=0; i<ncells; i++){
            particles[counter+0]->r = {.0+i};
            particles[counter+1]->r = {.5+i};
            counter += 2;
        }
    }
    else if(ndim == 2){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                particles[counter+0]->r = {.0+i, .0+j};
                particles[counter+1]->r = {.0+i, .5+j};
                particles[counter+2]->r = {.5+i, .0+j};
                counter += 4;
            }
        }
    }
    else if(ndim == 3){
        for(int i=0; i<ncells; i++){
            for(int j=0; j<ncells; j++){
                for(int k=0; k<ncells; k++){
                    particles[counter+0]->r = {.0+i, .0+j, .0+k};
                    particles[counter+1]->r = {.0+i, .5+j, .5+k};
                    particles[counter+2]->r = {.5+i, .0+j, .5+k};
                    particles[counter+3]->r = {.5+i, .5+j, .0+k};
                    counter += 4;
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
    for(Particle *particle : particles){
        particle->r *= scale;
    }

    return particles;
}

std::vector<class Particle *> from_xyz(std::string filename)
{
    std::vector<class Particle *> particles;  // = read_xyz(filename);
    return particles;
}
