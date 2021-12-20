#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

#include "io.h"
#include "box.h"
#include "dump.h"
#include "thermo.h"
#include "system.h"
#include "particle.h"
//#include "moves/moves.h"
#include "boundary/stillinger.h"
//#include "velocity/zero.h"


/* -------------------------------------------------------
   Box constructor 
---------------------------------------------------------- */

Box::Box(System* system_in)
{
    system = system_in;

    poteng = 0.;
    npar = 0;

    boundary = new Stillinger(this);
    //velocity = new Zero();

    // set default outputs
    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 1, "", outputs);
}


/* --------------------------------------------------------
   Set box boundaries
----------------------------------------------------------- */

void Box::set_boundary(class Boundary* boundary_in)
{
    boundary = boundary_in;
}


/* --------------------------------------------------
   Add a single particle from a particle object
----------------------------------------------------- */

void Box::add_particle(Particle* particle)
{
    npar ++;
    system->ndim = particle->r.size();
    particles.push_back(particle);
}


/* --------------------------------------------------
   Add a single particle given a label 'label' and
   initial position 'r'
----------------------------------------------------- */
   
void Box::add_particle(const std::string label, const std::valarray<double> r)
{
    npar ++;
    system->ndim = r.size();
    Particle *particle = new Particle(label, r);
    particles.push_back(particle);
}


/* --------------------------------------------------
   Add a set of particles, stored in a vector of
   particle objects 'particles_in'.
----------------------------------------------------- */

void Box::add_particles(std::vector<Particle *> particles_in)
{
    npar += particles_in.size();
    system->ndim = particles_in[0]->r.size();
    particles.insert(particles.end(), particles_in.begin(), particles_in.end());
}


/* ------------------------------------------------------
   Dump snapshot of system using the "write_xyz"-function. 
--------------------------------------------------------- */
   
void Box::snapshot(const std::string filename)
{
    std::vector<std::string> outputs = {"xyz"};
    Dump* tmp_dump = new Dump(this, 1, filename, outputs);
    tmp_dump->print_frame(0);
    delete tmp_dump;
}

/* -----------------------------------------------------
   Specify dump output
-------------------------------------------------------- */

void Box::set_dump(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    dump = new Dump(this, freq, filename, outputs);
}


/* -----------------------------------------------------
   Specify thermo output
-------------------------------------------------------- */

void Box::set_thermo(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    thermo = new Thermo(this, freq, filename, outputs);
}


/* -------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
---------------------------------------------------------------- */

std::vector<int> Box::build_neigh_list(const int i, const double rsq)
{
    double rijsq;
    std::valarray<double> ri = particles[i]->r;
    std::vector<int> neigh_list;
    for(int j=0; j<i; j++){
        rijsq = std::pow(particles[j]->r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for(int j=i+1; j<npar; j++){
        rijsq = std::pow(particles[j]->r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ----------------------------------------------------------
   Write number of times each system size has occured to
   file 'filename'
------------------------------------------------------------- */

void Box::write_nsystemsize(std::string filename)
{
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Gather(nsystemsize.data(), nsystemsize.size(), MPI_INT,  
    write_vector(nsystemsize, filename, "\n");
}
