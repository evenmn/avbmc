/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-09)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   This file hosts the Box class, and controls construction, manipulation and
   destruction of system boxes. Every box is associated with a boundary style
   and forcefield. In addition, particles are always affiliated with a box.
   Thermo and dump outputs are also associated with a box.
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

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
---------------------------------------------------------------------------- */

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
        std::cout << "memory_intensity has to be 1, 2 or 3! Aborting."
                  << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Copy constructor, needed to fulfill the rule of five needed for the class
   to be exception-safe.
---------------------------------------------------------------------------- */

Box::Box(const Box &box) : size_histogram(box.size_histogram),
    npartype(box.npartype), particles(box.particles),
    constraints(box.constraints)
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
---------------------------------------------------------------------------- */

void Box::swap(Box &other)
{
    unsigned int npar_tmp, step_tmp, ntype_tmp;
    unsigned int nmove_tmp, box_id_tmp, nconstraint_tmp;
    bool initialized_tmp, store_energy_tmp, store_distance_tmp;
    double time_tmp, poteng_tmp;
    std::vector<unsigned int> size_histogram_tmp, npartype_tmp;
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

    size_histogram_tmp = size_histogram;
    size_histogram = other.size_histogram;
    other.size_histogram = size_histogram_tmp;
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
---------------------------------------------------------------------------- */

void Box::set_forcefield(ForceField* forcefield_in)
{
    forcefield = forcefield_in;
    initialized = true;
    //if (store_energies) {
    //    forcefield->poteng_vec.resize(
}


/* ----------------------------------------------------------------------------
   Set box boundaries
---------------------------------------------------------------------------- */

void Box::set_boundary(Boundary* boundary_in)
{
    if (!boundary_allocated_externally) {
        delete boundary;
    }
    boundary = boundary_in;
    boundary_allocated_externally = true;
}


/* ----------------------------------------------------------------------------
   Add a single particle from a particle object. Particle is always appended
   to the particles-vector, with particle-ID npar-1. When adding a new
   particle, there are a few things that need to be done: 
     1. assign correct type based on the label
     2. add particle index to list of indices of particular type
     3. add particle to list of particles
     4. check that particle is located inside box, move it otherwise
     5. bump up number of particles and number of particles of particular type
     6. expand distance matrix and neighbor lists
     7. (expand size of energy and force matrices)
---------------------------------------------------------------------------- */

void Box::add_particle(Particle particle)
{
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before "
                  << "adding particles!" << std::endl;
        exit(0);
    }
    particle.type = forcefield->label2type[particle.label];
    typeidx[particle.type].push_back(npar);
    system->ndim = particle.r.size();
    particles.push_back(particle);
    boundary->correct_position(npar);
    npar ++;
    npartype[particle.type] ++;
    if (store_distance) {
        distance_manager->update_insert(npar-1);
    }
    //if (store_energy) {
    //    forcefield->update_insert(npar-1);
    //}
}


/* ----------------------------------------------------------------------------
   Add a single particle given a label 'label' and initial position 'r'
---------------------------------------------------------------------------- */

void Box::add_particle(const std::string &label, const std::valarray<double> r)
{
    add_particle(Particle(label, r));
}


/* ----------------------------------------------------------------------------
   Add a set of particles, stored in a vector of particle objects
   'particles_in'.
---------------------------------------------------------------------------- */

void Box::add_particles(std::vector<Particle> particles_in)
{
    for (Particle particle : particles_in) {
        add_particle(particle);
    }
}


/* ----------------------------------------------------------------------------
   Add a set of particles, stored in a vector of particle objects
   'particles_in'.
---------------------------------------------------------------------------- */

void Box::add_particles(const std::string &label,
    std::vector<std::valarray<double> > positions_in)
{
    for (std::valarray<double> position : positions_in) {
        add_particle(label, position);
    }
}


/* ----------------------------------------------------------------------------
   Initialize particles from an xyz-file, 'filename'
---------------------------------------------------------------------------- */

void Box::read_particles(const std::string &filename)
{
    std::vector<Particle> particles = read_xyz(filename);
    for (Particle particle : particles) {
        add_particle(particle);
    }
}


/* ----------------------------------------------------------------------------
   Remove a single particle by index. When a particle is removed, there are
   a few things that have to be updated:
     1. reduce number of particles and number of particles of particular 
        type by one
     2. remove particle index from list of indices of particular type
     3. reduce size of distance matrix and neighbor lists
     4. (reduce size of energy and force matrices)
---------------------------------------------------------------------------- */

void Box::_rm_typeidx(unsigned int i, unsigned int type)
{
    std::vector<unsigned int>::iterator position;

    position = std::find(typeidx[type].begin(), typeidx[type].end(), i);
    if (position != typeidx[type].end()) {
        typeidx[type].erase(position);
    }
    for (std::vector<unsigned int> &idxs : typeidx) {
        for (unsigned int &j : idxs) {
            if (j > i) {
                j--;
            }
        }
    }
}


void Box::rm_particle(unsigned int i)
{
    assert (i < npar);
    _rm_typeidx(i, particles[i].type);
    npar --;
    npartype[particles[i].type] --;
    if (store_distance) {
        distance_manager->update_remove(i);
    }
    //if (store_energy) {
    //    forcefield->update_remove(i);
    //}
    particles.erase(particles.begin() + i);
}


void Box::clear_particles()
{
    /*
    typeidx.clear();
    for (auto &type : npartype) {
        type.clear();
    }
    particles.clear();
    npar = 0;
    */
    for (unsigned int i=npar; i--;) {
        rm_particle(i);
    }
}


/* ----------------------------------------------------------------------------
   Add box constraint
---------------------------------------------------------------------------- */

void Box::add_constraint(Constraint* constraint)
{
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized "
                  << "before adding constraints!" << std::endl;
        exit(0);
    }
    std::cout << "Warning: All constraints have to be added before particles "
              << "are being added" << std::endl;
    nconstraint ++;
    constraints.push_back(constraint);
    constraint_allocated_in_system.push_back(false);
}


/* ----------------------------------------------------------------------------
   Remove box constraint
---------------------------------------------------------------------------- */

void Box::rm_constraint(unsigned int idx)
{
    if (constraint_allocated_in_system[idx]) {
        delete constraints[idx];
    }
    constraints.erase(constraints.begin() + idx);
    constraint_allocated_in_system.erase(
        constraint_allocated_in_system.begin() + idx
    );
    nconstraint --;
}


/* ----------------------------------------------------------------------------
   Dump snapshot of system using the "write_xyz"-function from io.cpp to file
   'filename', which is marked with box-ID and rank-ID if 'mark_file' is true
---------------------------------------------------------------------------- */
   
void Box::snapshot(std::string filename) //, const bool mark_file)
{
    std::vector<std::string> outputs = {"xyz"};
    Dump* dump_tmp = new Dump(this, 1, filename, outputs);
    dump_tmp->print_frame(0);
    delete dump_tmp;
}


/* ----------------------------------------------------------------------------
   Specify dump output
---------------------------------------------------------------------------- */

void Box::set_dump(const int freq, std::string filename, 
                   std::vector<std::string> outputs) //, const bool mark_file)
{
    if (!dump_allocated_externally) {
        delete dump;
    }
    dump = new Dump(this, freq, filename, outputs);
    dump_allocated_externally = true;
}


/* ----------------------------------------------------------------------------
   Specify thermo output
---------------------------------------------------------------------------- */

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
---------------------------------------------------------------------------- */
/*
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
*/

/* ----------------------------------------------------------------------------
   Write number of times each system size has occured to
   file 'filename'
---------------------------------------------------------------------------- */

void Box::write_size_histogram(const std::string &filename)
{
    unsigned int nsize;

    nsize = size_histogram.size();
    write_array(size_histogram.data(), nsize, filename, "\n");
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void Box::update_size_histogram()
{
    if (npar + 1 > size_histogram.size()) {
        size_histogram.resize(npar + 1);
    }
    size_histogram[npar] ++;
}


/* ----------------------------------------------------------------------------
   Box destructor, deleting boundary, distance_manager, thermo and dump
   pointers
---------------------------------------------------------------------------- */

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
