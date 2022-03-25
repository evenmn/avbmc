#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "particle.h"

class Box
{
public:
    Box(class System *, int = 2); 

    // methods
    void set_forcefield(class ForceField*);
    void set_boundary(class Boundary *);
    void add_particle(Particle);
    void add_particle(std::string, std::valarray<double>);
    void add_particles(std::vector<Particle>);
    void add_particles(std::string, std::vector<std::valarray<double> >);
    void add_constraint(class Constraint *);

    std::string file_marking();
    void snapshot(std::string, bool = true);
    void set_dump(int, std::string, std::vector<std::string>, bool = true);
    void set_thermo(int, std::string, std::vector<std::string>, bool = true);
    std::vector<int> build_neigh_list(int, double);
    std::vector<int> build_neigh_list(int, double**);
    void write_nsystemsize(std::string);
    ~Box();

    // variables
    class System* system = nullptr;
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;
    class ForceField* forcefield = nullptr;
    //class Velocity* velocity = nullptr;
    class DistanceManager* distance_manager = nullptr;

    bool initialized, store_energy, store_distance;
    unsigned int npar, step, ntype, nmove, box_id, nconstraint;
    double poteng, time;

    std::vector<int> nsystemsize, npartype;
    std::vector<Particle> particles;
    std::vector<class Constraint *> constraints;
};
