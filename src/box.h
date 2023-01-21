/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "particle.h"

class Box
{
public:
    Box(class System *); //, int = 2); 
    Box(const Box &);  // copy constructor
    Box& operator=(const Box &other) // overloading assignment operator
    {
        Box tmp(other); // calling copy constructor
        swap(tmp);
        return *this;
    };
    void swap(Box &other);  // swap two objects

    // methods
    void set_forcefield(class ForceField*);
    void set_boundary(class Boundary *);
    //void set_integrator(class Integrator *);
    void add_particle(Particle);
    void add_particle(const std::string &, const std::valarray<double> &);
    void add_particles(std::vector<Particle> &);
    void add_particles(const std::string &, std::vector<std::valarray<double> > &);
    void rm_particle(const unsigned int &);
    void clear_particles();
    void read_particles(const std::string &);
    void add_constraint(class Constraint *);
    void _rm_typeidx(const unsigned int &, const unsigned int &);
    void rm_constraint(unsigned int);

    std::string file_marking();
    void snapshot(std::string); //, bool = true);
    void set_dump(int, std::string, std::vector<std::string>); //, bool = true);
    void set_thermo(int, std::string, std::vector<std::string>); //, bool = true);
    //std::vector<unsigned int> build_neigh_list(int, double);
    //std::vector<unsigned int> build_neigh_list(int, double**);
    void update_size_histogram();
    void write_size_histogram(const std::string &);
    ~Box();

    // variables
    class System* system = nullptr;
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;
    class ForceField* forcefield = nullptr;
    //class Integrator* integrator = nullptr;
    //class Velocity* velocity = nullptr;
    class DistanceManager* distance_manager = nullptr;

    bool initialized, store_energy, store_distance, boundary_allocated_externally;
    bool dump_allocated_externally, thermo_allocated_externally;
    bool boundary_allocated_in_system, forcefield_allocated_in_system, box_allocated_in_system;
    unsigned int npar, step, ntype, nmove, box_id, nconstraint;
    double poteng, time;

    std::vector<unsigned int> size_histogram, npartype;  // npartype is used by constraints
    std::vector<std::vector<unsigned int> > typeidx;     // indices of each type
    std::vector<Particle> particles;
    std::vector<class Constraint *> constraints;
    std::vector<bool> constraint_allocated_in_system;
};
