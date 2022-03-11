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
#include "distance_manager.h"
#include "forcefield/forcefield.h"
#include "boundary/open.h"
//#include "velocity/zero.h"


/* ----------------------------------------------------------------------------
   Box constructor. The argument 'memory_intensity' specifies how memeory
   intensive the simulation should be, with a memory-cpu-time tradeoff. The
   options are:

       1: Storing necessary neighbor lists only
       2: Storing distances and relative coordinates between particles
       3: Storing distances, relative coordinates and energy contributions
          of each particle
------------------------------------------------------------------------------- */

Box::Box(System* system_in, const int memory_intensity)
{
    system = system_in;

    poteng = 0.;
    initialized = false;
    npar = ntype = nmove = step = nconstraint = 0;

    //velocity = new Zero();
    distance_manager = new DistanceManager(this);

    // set default outputs
    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 0, "", outputs);

    // memory intensitivity
    if (memory_intensity==1) {
        store_distance = false;
        store_energy = false;
    }
    else if (memory_intensity==2) {
        store_distance = true;
        store_energy = false;
    }
    else if (memory_intensity==3) {
        store_distance = true;
        store_energy = true;
    }
    else {
        std::cout << "memory_intensity has to be 1, 2 or 3! Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Overwrite default forcefield object
------------------------------------------------------------------------------- */

void Box::set_forcefield(ForceField* forcefield_in)
{
    forcefield = forcefield_in;
    initialized = true;
    //if (store_energies) {
    //    forcefield->poteng_vec.resize(
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
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before adding particles!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
    npar ++;
    particle.type = forcefield->label2type.at(particle.label);
    npartype.resize(forcefield->ntype);
    npartype[particle.type] ++;
    system->ndim = particle.r.size();
    particles.push_back(particle);
}


/* ----------------------------------------------------------------------------
   Add a single particle given a label 'label' and initial position 'r'
------------------------------------------------------------------------------- */

void Box::add_particle(const std::string label, const std::valarray<double> r)
{
    add_particle(Particle(label, r));
}


/* ----------------------------------------------------------------------------
   Add a set of particles, stored in a vector of particle objects
   'particles_in'.
------------------------------------------------------------------------------- */

void Box::add_particles(std::vector<Particle> particles_in)
{
    for (Particle particle : particles_in) {
        add_particle(particle);
    }
}


/* ----------------------------------------------------------------------------
   Add box constraint
------------------------------------------------------------------------------- */

void Box::add_constraint(class Constraint* constraint)
{
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before adding constraints!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
    nconstraint ++;
    constraints.push_back(constraint);
}


/* ----------------------------------------------------------------------------
   Returning a string with information about rank-ID and box-ID if there are
   more than one rank/box. The string can then be added to a filename to mark
   the file.
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
        //std::valarray<double> rii = particles[j].r - particles[i].r;
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}

std::vector<int> Box::build_neigh_list(const int i, double **rsq)
{
    unsigned int typei, typej, equaltypecount;
    double rijsq;
    typei = particles[i].type;
    std::valarray<double> ri = particles[i].r;
    std::vector<int> neigh_list;
    for(unsigned int j=0; j<i; j++){
        typej = particles[j].type;
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq[typei][typej]){
            neigh_list.push_back(j);
            if (typej == typei) equaltypecount++;
        }
    }
    for(unsigned int j=i+1; j<npar; j++){
        typej = particles[j].type;
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq[typei][typej]){
            neigh_list.push_back(j);
            if (typej == typei) equaltypecount++;
        }
    }
    if (particles[i].label == "O" && npar > 6 && equaltypecount < 2) {
        neigh_list.clear();
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
    delete distance_manager;
}
