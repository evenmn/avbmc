#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

#include <mpi.h>

#include "particle.h"
#include "io.h"
#include "box.h"
#include "dump.h"
#include "thermo.h"
#include "system.h"
#include "forcefield/forcefield.h"
#include "boundary/stillinger.h"
//#include "velocity/zero.h"


/* ----------------------------------------------------------------------------
   Box constructor 
------------------------------------------------------------------------------- */

Box::Box(System* system_in)
{
    system = system_in;

    poteng = 0.;
    npar = ntype = nmove = step = 0;

    //velocity = new Zero();

    // set default outputs
    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 0, "", outputs);
}


/* ----------------------------------------------------------------------------
   Set box boundaries
------------------------------------------------------------------------------- */

void Box::set_boundary(class Boundary* boundary_in)
{
    boundary = boundary_in;
}


/* ----------------------------------------------------------------------------
   Add a single particle from a particle object
------------------------------------------------------------------------------- */

void Box::add_particle(Particle particle)
{
    if (!system->initialized) {
        std::cout << "Forcefield needs to be initialized before adding particles!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
    npar ++;
    particle.type = system->forcefield->label2type.at(particle.label);
    system->ndim = particle.r.size();
    particles.push_back(particle);
}


/* ----------------------------------------------------------------------------
   Add a single particle given a label 'label' and
   initial position 'r'
------------------------------------------------------------------------------- */

void Box::add_particle(const std::string label, const std::valarray<double> r)
{
    add_particle(Particle(label, r));
}


/* ----------------------------------------------------------------------------
   Add a set of particles, stored in a vector of
   particle objects 'particles_in'.
------------------------------------------------------------------------------- */

void Box::add_particles(std::vector<Particle> particles_in)
{
    for (Particle particle : particles_in) {
        add_particle(particle);
    }
}


/* ----------------------------------------------------------------------------
   Returning a string with information about rank-ID and box-ID if there are
   more than one rank/box. The string can then be added to a filename to 
   mark the file.
------------------------------------------------------------------------------- */

std::string Box::file_marking()
{
    std::string str;
    if (system->nbox > 1) {
        str = "BOX" + std::to_string(box_id) + "_";
    }
    if (system->nprocess > 1) {
        str += "RANK" + std::to_string(system->rank) + "_";
    }
    return str;
}

/* ----------------------------------------------------------------------------
   Dump snapshot of system using the "write_xyz"-function from io.cpp to file
   'filename', which is marked with box-ID and rank-ID if 'mark_file' is true
------------------------------------------------------------------------------- */
   
void Box::snapshot(std::string filename, const bool mark_file)
{
    if (mark_file) {
        filename = file_marking() + filename;
    }
    std::vector<std::string> outputs = {"xyz"};
    Dump* tmp_dump = new Dump(this, 1, filename, outputs);
    tmp_dump->print_frame(0);
    delete tmp_dump;
}

/* ----------------------------------------------------------------------------
   Specify dump output
------------------------------------------------------------------------------- */

void Box::set_dump(const int freq, std::string filename, 
                   std::vector<std::string> outputs, const bool mark_file)
{
    if (mark_file) {
        filename = file_marking() + filename;
    }
    dump = new Dump(this, freq, filename, outputs);
}


/* ----------------------------------------------------------------------------
   Specify thermo output
------------------------------------------------------------------------------- */

void Box::set_thermo(const int freq, std::string filename,
                     std::vector<std::string> outputs, const bool mark_file)
{
    if (mark_file) {
        filename = file_marking() + filename;
    }
    thermo = new Thermo(this, freq, filename, outputs);
}


/* ----------------------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
------------------------------------------------------------------------------- */

double normsq(std::valarray<double> array)
{
    double sumsq = 0.;
    for(double element : array){
        sumsq += element * element;
    }
    return sumsq;
}

std::vector<int> Box::build_neigh_list(const int i, const double rsq)
{
    double rijsq;
    std::valarray<double> ri = particles[i].r;
    std::vector<int> neigh_list;
    for(unsigned int j=0; j<i; j++){
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for(unsigned int j=i+1; j<npar; j++){
        std::valarray<double> rii = particles[j].r - particles[i].r;
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ----------------------------------------------------------------------------
   Write number of times each system size has occured to
   file 'filename'
------------------------------------------------------------------------------- */

void Box::write_nsystemsize(std::string filename)
{
    int maxsize;
    MPI_Barrier(MPI_COMM_WORLD);
    int size = nsystemsize.size();
    MPI_Reduce(&size, &maxsize, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nsystemsize.resize(maxsize);
    int* nsystemsizetot = new int[maxsize];
    MPI_Reduce(nsystemsize.data(), nsystemsizetot, maxsize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (system->rank == 0)
    {
        write_array(nsystemsizetot, maxsize, filename, "\n");
    }
    delete[] nsystemsizetot;
}


/* ----------------------------------------------------------------------------
   Box destructor, deleting thermo and dump pointers
------------------------------------------------------------------------------- */
Box::~Box()
{
    delete dump;
    delete thermo;
}
