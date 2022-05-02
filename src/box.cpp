#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

//#include <mpi.h>

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

Box::Box(System* system_in) //, int memory_intensity)
{
    system = system_in;

    time = poteng = 0.;
    initialized = boundary_allocated_externally = dump_allocated_externally =
        thermo_allocated_externally = box_allocated_in_system = false;
    boundary_allocated_in_system = forcefield_allocated_in_system = false;
    npar = ntype = nmove = step = nconstraint = box_id = 0;

    // set default values
    boundary = new Open(this);
    //velocity = new Zero();
    distance_manager = new DistanceManager(this);

    // set default outputs
    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 0, "", outputs);

    // memory intensitivity
    int memory_intensity = 2;
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
   Copy constructor, needed to fulfill the rule of five needed for the class
   to be exception-safe.
------------------------------------------------------------------------------- */

Box::Box(const Box &box) : nsystemsize(box.nsystemsize), npartype(box.npartype),
                           particles(box.particles), constraints(box.constraints)
{
    //system = new System(box.system);
    //dump = new Dump();
    //thermo = new Thermo();
    //boundary = new Boundary();
    //forcefield = new Forcefield();
    //integrator = new integrator();
    //velocity = new velocity();
    //distance_manager = new DistanceManager();

    //*system = *box.system;
    //*dump = *box.dump;
    //*thermo = *box.thermo;
    //*boundary = *box.boundary;
    //*forcefield = *box.forcefield;
    //*integrator = *box.integrator;
    //*velocity = *box.velocity;
    *distance_manager = *box.distance_manager;

    initialized = box.initialized;
    store_energy = box.store_energy;
    store_distance = box.store_distance;
    npar = box.npar;
    step = box.step;
    ntype = box.ntype;
    nmove = box.nmove;
    box_id = box.box_id;
    nconstraint = box.nconstraint;
    poteng = box.poteng;
    time = box.time;
}


/* ----------------------------------------------------------------------------
   Swap the ownership of the internals
------------------------------------------------------------------------------- */

void Box::swap(Box &other)
{
    unsigned int npar_tmp, step_tmp, ntype_tmp, nmove_tmp, box_id_tmp, nconstraint_tmp;
    bool initialized_tmp, store_energy_tmp, store_distance_tmp;
    double time_tmp, poteng_tmp;
    std::vector<unsigned int> nsystemsize_tmp, npartype_tmp;
    std::vector<Particle> particles_tmp;
    std::vector<Constraint *> constraints_tmp;

    System *system_tmp = system;
    system = other.system;
    other.system = system_tmp;
    Dump *dump_tmp = dump;
    dump = other.dump;
    other.dump = dump_tmp;
    Thermo *thermo_tmp = thermo;
    thermo = other.thermo;
    other.thermo = thermo_tmp;
    Boundary *boundary_tmp = boundary;
    boundary = other.boundary;
    other.boundary = boundary_tmp;
    ForceField *forcefield_tmp = forcefield;
    forcefield = other.forcefield;
    other.forcefield = forcefield_tmp;
    //integrator = box.integrator;
    //velocity = box.velocity;
    DistanceManager *distance_manager_tmp = distance_manager;
    distance_manager = other.distance_manager;
    other.distance_manager = distance_manager_tmp;

    initialized_tmp = initialized;
    initialized = other.initialized;
    other.initialized = initialized_tmp;
    store_energy_tmp = store_energy;
    store_energy = other.store_energy;
    other.store_energy = store_energy_tmp;
    store_distance_tmp = store_distance;
    store_distance = other.store_distance;
    other.store_distance = store_distance_tmp;

    npar_tmp = npar;
    npar = other.npar;
    other.npar = npar_tmp;
    step_tmp = step;
    step = other.step;
    other.step = step_tmp;
    ntype_tmp = ntype;
    ntype = other.ntype;
    other.ntype = ntype_tmp;
    nmove_tmp = nmove;
    nmove = other.nmove;
    other.nmove = nmove_tmp;
    box_id_tmp = box_id;
    box_id = other.box_id;
    other.box_id = box_id_tmp;
    nconstraint_tmp = nconstraint;
    nconstraint = other.nconstraint;
    other.nconstraint = nconstraint_tmp;

    poteng_tmp = poteng;
    poteng = other.poteng;
    other.poteng = poteng_tmp;
    time_tmp = time;
    time = other.time;
    other.time = time_tmp;

    nsystemsize_tmp = nsystemsize;
    nsystemsize = other.nsystemsize;
    other.nsystemsize = nsystemsize_tmp;
    npartype_tmp = npartype;
    npartype = other.npartype;
    other.npartype = npartype_tmp;
    particles_tmp = particles;
    particles = other.particles;
    other.particles = particles_tmp;
    constraints_tmp = constraints;
    constraints = other.constraints;
    other.constraints = constraints_tmp;
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

void Box::set_boundary(Boundary* boundary_in)
{
    if (!boundary_allocated_externally) {
        delete boundary;
    }
    boundary = boundary_in;
    boundary_allocated_externally = true;
}


/* ----------------------------------------------------------------------------
   Add a single particle from a particle object
------------------------------------------------------------------------------- */

void Box::add_particle(Particle particle)
{
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before adding particles!" << std::endl;
        //MPI_Abort(MPI_COMM_WORLD, 143);
        exit(0);
    }
    npar ++;
    particle.type = forcefield->label2type[particle.label];
    npartype.resize(forcefield->ntype);
    npartype[particle.type] ++;
    system->ndim = particle.r.size();
    particles.push_back(particle);
}


/* ----------------------------------------------------------------------------
   Add a single particle given a label 'label' and initial position 'r'
------------------------------------------------------------------------------- */

void Box::add_particle(const std::string &label, const std::valarray<double> r)
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
   Add a set of particles, stored in a vector of particle objects
   'particles_in'.
------------------------------------------------------------------------------- */

void Box::add_particles(const std::string &label, std::vector<std::valarray<double> > positions_in)
{
    for (std::valarray<double> position : positions_in) {
        add_particle(label, position);
    }
}


/* ----------------------------------------------------------------------------
   Initialize particles from an xyz-file, 'filename'
------------------------------------------------------------------------------- */

void Box::read_particles(const std::string &filename)
{
    std::vector<Particle> particles = read_xyz(filename);
    for (Particle particle : particles) {
        add_particle(particle);
    }
}


/* ----------------------------------------------------------------------------
   Add box constraint
------------------------------------------------------------------------------- */

void Box::add_constraint(Constraint* constraint)
{
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before adding constraints!" << std::endl;
        //MPI_Abort(MPI_COMM_WORLD, 143);
        exit(0);
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
   
void Box::snapshot(std::string filename) //, const bool mark_file)
{
    //if (mark_file) {
    //    filename = file_marking() + filename;
    //}
    std::vector<std::string> outputs = {"xyz"};
    Dump* dump_tmp = new Dump(this, 1, filename, outputs);
    dump_tmp->print_frame(0);
    delete dump_tmp;
}

/* ----------------------------------------------------------------------------
   Specify dump output
------------------------------------------------------------------------------- */

void Box::set_dump(const int freq, std::string filename, 
                   std::vector<std::string> outputs) //, const bool mark_file)
{
    //if (mark_file) {
    //    filename = file_marking() + filename;
    //}
    if (!dump_allocated_externally) {
        delete dump;
    }
    dump = new Dump(this, freq, filename, outputs);
    dump_allocated_externally = true;
}


/* ----------------------------------------------------------------------------
   Specify thermo output
------------------------------------------------------------------------------- */

void Box::set_thermo(const int freq, std::string filename,
                     std::vector<std::string> outputs) //, const bool mark_file)
{
    //if (mark_file) {
    //    filename = file_marking() + filename;
    //}
    if (!thermo_allocated_externally) {
        delete thermo;
    }
    thermo = new Thermo(this, freq, filename, outputs);
    thermo_allocated_externally = true;
}


/* ----------------------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
------------------------------------------------------------------------------- */

double normsq(std::valarray<double> array)
{
    double sumsq;

    sumsq = 0.;
    for(double element : array){
        sumsq += element * element;
    }
    return sumsq;
}

std::vector<unsigned int> Box::build_neigh_list(const int i, const double rsq)
{
    double rijsq;
    std::valarray<double> ri = particles[i].r;
    std::vector<unsigned int> neigh_list;
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

std::vector<unsigned int> Box::build_neigh_list(const int i, double **rsq)
{
    unsigned int typei, typej, j, equaltypecount;
    double rijsq;
    typei = particles[i].type;
    std::valarray<double> ri = particles[i].r;
    std::vector<unsigned int> neigh_list;

    equaltypecount = 0;
    for(j=0; j<i; j++){
        typej = particles[j].type;
        rijsq = normsq(particles[j].r - ri);
        if(rijsq < rsq[typei][typej]){
            neigh_list.push_back(j);
            if (typej == typei) equaltypecount++;
        }
    }
    for(j=i+1; j<npar; j++){
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

void Box::write_nsystemsize(const std::string &filename)
{
    //int maxsize;
    //MPI_Barrier(MPI_COMM_WORLD);
    unsigned int size = nsystemsize.size();
    //MPI_Reduce(&size, &maxsize, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&maxsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //nsystemsize.resize(maxsize);
    //int* nsystemsizetot = new int[maxsize];
    //MPI_Reduce(nsystemsize.data(), nsystemsizetot, maxsize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (system->rank == 0)
    //{
    //    write_array(nsystemsizetot, maxsize, filename, "\n");
    //}
    //delete[] nsystemsizetot;
    write_array(nsystemsize.data(), size, filename, "\n");
}


/* ----------------------------------------------------------------------------
   Box destructor, deleting thermo and dump pointers
------------------------------------------------------------------------------- */
Box::~Box()
{
    if (!dump_allocated_externally) {
        delete dump;
    }
    if (!thermo_allocated_externally) {
        delete thermo;
    }
    if (!boundary_allocated_externally) {
        delete boundary;
    }
    delete distance_manager;
}
