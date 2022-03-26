#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

#include "moves.h"
#include "../particle.h"
#include "../system.h"
#include "../rng/rng.h"


/* ----------------------------------------------------------------------------
   Moves base class constructor
------------------------------------------------------------------------------- */

Moves::Moves(System* system_in)
{
    system = system_in;
    rng = system->rng;
}


/* ----------------------------------------------------------------------------
   Rotate molecule around the center of mass atom by a random angle. 
------------------------------------------------------------------------------- */

std::vector<Particle> Moves::rotate_molecule(std::vector<Particle> particles)
{
    if(system->ndim == 1){
        // cannot rotate in 1D
    }
    else if(system->ndim == 2){
        double angle = 2 * pi * rng->next_double();
        for(Particle &particle : particles){
            std::valarray<double> r = particle.r;
            particle.r = {
                r[0] * std::cos(angle) - r[1] * std::sin(angle),
                r[0] * std::sin(angle) + r[1] * std::cos(angle)
            };
        }
    }
    else{
        double anglea = 2 * pi * rng->next_double();
        double cosa = std::cos(anglea);
        double sina = std::sin(anglea);

        double angleb = 2 * pi * rng->next_double();
        double cosb = std::cos(angleb);
        double sinb = std::sin(angleb);

        double anglec = 2 * pi * rng->next_double();
        double cosc = std::cos(anglec);
        double sinc = std::sin(anglec);

        double Axx = cosa * cosb;
        double Axy = cosa * sinb * sinc - sina * cosc;
        double Axz = cosa * sinb * cosc + sina * sinc;

        double Ayx = sina * cosb;
        double Ayy = sina * sinb * sinc + cosa * cosc;
        double Ayz = sina * sinb * cosc - cosa * sinc;

        double Azx = -sinb;
        double Azy = cosb * sinc;
        double Azz = cosb * cosc;

        for (Particle &particle : particles) {
            std::valarray<double> r = particle.r;
            particle.r[0] = Axx*r[0] + Axy*r[1] + Axz*r[2];
            particle.r[1] = Ayx*r[0] + Ayy*r[1] + Ayz*r[2];
            particle.r[2] = Azx*r[0] + Azy*r[1] + Azz*r[2];
        }
    }
    return particles;
}


/* ----------------------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
------------------------------------------------------------------------------- */

double Moves::norm(std::valarray<double> array)
{
    double normsq;
    unsigned int i;
    normsq = 0.;
    for (i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


/* ----------------------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor distance squared
   'rsq'
------------------------------------------------------------------------------- */

std::vector<int> Moves::build_neigh_list(std::vector<Particle> particles,
                                         const int i, const double rsq)
{
    double rijsq;
    unsigned int npar, j;
    
    npar = particles.size();
    std::valarray<double> ri = particles[i].r;
    std::vector<int> neigh_list;
    for (j=0; j<i; j++) {
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        //std::cout << "rij " << rsq << " " << rijsq << " " << norm(particles[j].r - ri) << std::endl;
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for (j=i+1; j<npar; j++) {
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        //std::cout << "rij " << rsq << " " << rijsq << " " << norm(particles[j].r - ri) << std::endl;
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


/* ----------------------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------------------- */

void Moves::check_neigh_recu(const int i, std::vector<Particle> molecule,
                            unsigned int elm_count, std::vector<int> &elm_idx,
                            std::vector<Particle> particles, double rc) {
    //std::cout << "1elm_idx.size() " << elm_idx.size() << std::endl;
    if (elm_idx.size() < molecule.size()) {  
        //std::cout << "2elm_idx.size() " << elm_idx.size() << std::endl;
        //std::cout << particles[i].type << " " << molecule[elm_count].type << std::endl;
        if (particles[i].type == molecule[elm_count].type) {
            //std::cout << "1molecule.size(): " << elm_idx.size() << " " << molecule.size() << std::endl;
            elm_idx.push_back(i);
            elm_count ++;
            //std::cout << "2particles.size() " << particles.size() << std::endl;
            std::vector<int> neigh_list = build_neigh_list(particles, i, rc*rc);
            //std::cout << "3neigh_list.size(): " << neigh_list.size() << std::endl;
            for (int neigh : neigh_list) {
                check_neigh_recu(neigh, molecule, elm_count, elm_idx, particles, rc);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
------------------------------------------------------------------------------- */

std::vector<int> Moves::detect_molecule(std::vector<Particle> particles,
                                        std::vector<Particle> molecule,
                                        bool &detected, double rc)
{
    unsigned int i, count;
    std::vector<int> elm_idx;

    /*
    std::cout << "--- molecule types: ";
    for (Particle particle : molecule) {
        std::cout << particle.type << " ";
    }
    std::cout << std::endl;
    
    std::cout << "--- particle types: ";
    for (Particle particle : particles) {
        std::cout << particle.type << " ";
    }
    std::cout << std::endl;
    */

    count = 0;
    while (count < particles.size() && !detected)
    {
        elm_idx.clear();
        i = rng->next_int(particles.size());     // pick initial particle
        check_neigh_recu(i, molecule, 0, elm_idx, particles, rc);
        if (elm_idx.size() == molecule.size()) {
            detected = true;
        }
        count ++;
    }
    if (!detected)
    {
        elm_idx.clear();
    }
    return elm_idx;
}


/* ----------------------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------------------- */

void Moves::check_neigh_recu(const int i, std::vector<Particle> molecule,
                            unsigned int elm_count, std::vector<int> &elm_idx,
                            std::vector<std::vector<int> > neigh_list) {
    if (elm_idx.size() < molecule.size()) {  
        //if (particles[i].label == molecule[elm_count].label) {
        elm_idx.push_back(i);
        elm_count ++;
        for (int neigh : neigh_list[i]) {
            check_neigh_recu(neigh, molecule, elm_count, elm_idx, neigh_list);
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
------------------------------------------------------------------------------- */

std::vector<int> Moves::detect_molecule(std::vector<std::vector<int> > neigh_list,
                                        std::vector<Particle> molecule,
                                        bool &detected)
{
    unsigned int i, count;
    std::vector<int> elm_idx;
    
    count = 0;
    while (count < neigh_list.size() && !detected)
    {
        elm_idx.clear();
        i = rng->next_int(neigh_list.size());     // pick initial particle
        check_neigh_recu(i, molecule, 0, elm_idx, neigh_list);
        if (elm_idx.size() == molecule.size()) {
            detected = true;
        }
        count ++;
    }
    if (!detected)
    {
        elm_idx.clear();
    }
    return elm_idx;
}


/* ----------------------------------------------------------------------------
   Check if particles match molecule type recursively.
------------------------------------------------------------------------------- */

void Moves::check_neigh_recu(const int i, std::vector<Particle> particles,
                            std::vector<Particle> molecule,
                            unsigned int elm_count, std::vector<int> &elm_idx,
                            std::vector<std::vector<int> > neigh_list) {
    if (elm_idx.size() < molecule.size()) {  
        //std::cout << "1molecule.size(): " << elm_idx.size() << " " << molecule.size() << std::endl;
        if (particles[i].type == molecule[elm_count].type) {
            //std::cout << "2molecule.size(): " << elm_idx.size() << " " << molecule.size() << std::endl;
            elm_idx.push_back(i);
            elm_count ++;
            //std::cout << "neigh_list.size(): " << neigh_list[i].size() << std::endl;
            for (int neigh : neigh_list[i]) {
                check_neigh_recu(neigh, particles, molecule, elm_count, elm_idx, neigh_list);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
------------------------------------------------------------------------------- */

std::vector<int> Moves::detect_molecule(std::vector<std::vector<int> > neigh_list,
                                        std::vector<Particle> particles, 
                                        std::vector<Particle> molecule,
                                        bool &detected)
{
    unsigned int i, count;
    std::vector<int> elm_idx;
    
    count = 0;
    while (count < neigh_list.size() && !detected)
    {
        elm_idx.clear();
        i = rng->next_int(neigh_list.size());     // pick initial particle
        check_neigh_recu(i, particles, molecule, 0, elm_idx, neigh_list);
        if (elm_idx.size() == molecule.size()) {
            detected = true;
        }
        count ++;
    }
    if (!detected)
    {
        elm_idx.clear();
    }
    return elm_idx;
}


