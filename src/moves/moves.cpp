/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  This file hosts the move base class. It has methods for finding insertion
  positions of particles, rotating and detecting molecules, among others.
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>
#include <cassert>

#include "moves.h"
#include "../box.h"
#include "../particle.h"
#include "../system.h"
#include "../rng/rng.h"
#include "../distance_manager.h"


/* ----------------------------------------------------------------------------
   Moves base class constructor
---------------------------------------------------------------------------- */

Moves::Moves(System* system_in)
{
    system = system_in;
    rng = system->rng;
    ndrawn = naccept = 0;
    cum_time = du = 0.;
}


/* ----------------------------------------------------------------------------
   Rotate molecule around the center of mass atom by a random angle. 
---------------------------------------------------------------------------- */

std::vector<Particle> Moves::rotate_molecule(std::vector<Particle> particles)
{
    if (system->ndim == 1) {
        // cannot rotate in 1D
    }
    else if (system->ndim == 2) {
        double angle = 2 * pi * rng->next_double();
        for(Particle &particle : particles){
            std::valarray<double> r = particle.r;
            particle.r = {
                r[0] * std::cos(angle) - r[1] * std::sin(angle),
                r[0] * std::sin(angle) + r[1] * std::cos(angle)
            };
        }
    }
    else {
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
---------------------------------------------------------------------------- */

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
   Find where to insert particle relative to target particle. There are two
   common techniques for sampling positions uniformly from a hypersphere:
    1. Sample hypercube and reject points that are not inside hypersphere
    2. Obtain positions in spherical coordinates and transform to Cartesian
   The former method is recommended for low dimensions (d<4), as less than
   50% of attempts will be rejected in average. Set spherical=false to use it.
---------------------------------------------------------------------------- */

std::valarray<double> Moves::insertion_position(bool spherical)
{
    std::valarray<double> dr(system->ndim);
    if (spherical) {
        double r, phi, cos_theta, pref;

        assert (system->ndim == 3);
        phi = 2 * pi * rng->next_double();
        r = std::pow(rng->next_double(), 1/3.) * r_above;
        cos_theta = 2 * rng->next_double() - 1.;
        pref = r * std::sqrt(1. - cos_theta * cos_theta);
        dr[0] = pref * std::cos(phi);
        dr[1] = pref * std::sin(phi);
        dr[2] = r * cos_theta;
    }
    else {
        double normsq;

        normsq = norm(dr);
        while (normsq > r_abovesq || normsq < r_belowsq) {
            for (double &d : dr) {
                d = r_above * (2 * rng->next_double() - 1);
            }
            normsq = norm(dr);
        }
    }
    return dr;
}


/* ----------------------------------------------------------------------------
   Detect target particle. Doing npar attempts of detecting particle. This
   could be done much more efficiently by constructing a look up table with
   particle indices of each type.
---------------------------------------------------------------------------- */
/*
unsigned int Moves::detect_target_particle(bool &detected)
{
    unsigned int i, count;

    detected = false;
    count = 0;
    while (!detected && count < box->npar) {
        i = rng->next_int(box->npar);
        if (box->particles[i].label == particle_label) {
            detected = true;
        }
        count ++;
    }
    return i;
}
*/

/* ----------------------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor distance squared
   'rsq'
---------------------------------------------------------------------------- */
/*
std::vector<unsigned int> Moves::build_neigh_list(std::vector<Particle> particles,
    const unsigned int i, const double rsq)
{
    double rijsq;
    unsigned int npar, j;
    
    npar = particles.size();
    std::valarray<double> ri = particles[i].r;
    std::vector<unsigned int> neigh_list;
    for (j=0; j<i; j++) {
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        if(rijsq < rsq){
            neigh_list.push_back(j);
        }
    }
    for (j=i+1; j<npar; j++) {
        rijsq = std::pow(particles[j].r - ri, 2).sum();
        if (rijsq < rsq) {
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}
*/

/* ----------------------------------------------------------------------------
   Check if particles match molecule type recursively.
---------------------------------------------------------------------------- */

void Moves::check_neigh_recu(const int i, std::vector<Particle> molecule,
    unsigned int elm_count, std::vector<unsigned int> &elm_idx,
    std::vector<Particle> particles, double rc, Box *box)
{
    if (elm_idx.size() < molecule.size()) {  
        if (particles[i].type == molecule[elm_count].type) {
            elm_idx.push_back(i);
            elm_count ++;
            std::vector<unsigned int> neigh_list = box->distance_manager->build_neigh_list(particles, i, rc*rc);
            for (int neigh : neigh_list) {
                check_neigh_recu(neigh, molecule, elm_count, elm_idx, particles, rc, box);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
---------------------------------------------------------------------------- */

std::vector<unsigned int> Moves::detect_molecule(std::vector<Particle> particles,
    std::vector<Particle> molecule, bool &detected, double rc, Box *box)
{
    unsigned int i, count;
    std::vector<unsigned int> elm_idx;

    count = 0;
    while (count < particles.size() && !detected)
    {
        elm_idx.clear();
        i = rng->next_int(particles.size());     // pick initial particle
        check_neigh_recu(i, molecule, 0, elm_idx, particles, rc, box);
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
---------------------------------------------------------------------------- */

void Moves::check_neigh_recu(const int i, std::vector<Particle> molecule,
    unsigned int elm_count, std::vector<unsigned int> &elm_idx,
    std::vector<std::vector<unsigned int> > neigh_list)
{
    if (elm_idx.size() < molecule.size()) {  
        //if (particles[i].label == molecule[elm_count].label) {
        elm_idx.push_back(i);
        elm_count ++;
        for (unsigned int neigh : neigh_list[i]) {
            check_neigh_recu(neigh, molecule, elm_count, elm_idx, neigh_list);
        }
    }
}


/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
---------------------------------------------------------------------------- */

std::vector<unsigned int> Moves::detect_molecule(
    std::vector<std::vector<unsigned int> > neigh_list,
    std::vector<Particle> molecule, bool &detected)
{
    unsigned int i, count;
    std::vector<unsigned int> elm_idx;
    
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
---------------------------------------------------------------------------- */
/*
void Moves::check_neigh_recu(const int i, std::vector<Particle> particles,
    std::vector<Particle> molecule, unsigned int elm_count,
    std::vector<unsigned int> &elm_idx,
    std::vector<std::vector<unsigned int> > neigh_list)
{
    if (elm_idx.size() < molecule.size()) {  
        if (particles[i].type == molecule[elm_count].type) {
            elm_idx.push_back(i);
            elm_count ++;
            for (unsigned int neigh : neigh_list[i]) {
                check_neigh_recu(neigh, particles, molecule, elm_count,
                    elm_idx, neigh_list);
            }
        }
    }
}
*/

/* ----------------------------------------------------------------------------
   Detect molecule of the same types as 'molecule' randomly by picking a random 
   atom among the elements and checking the neighbor list.
   Returning a list of atom ids if molecule is detected
---------------------------------------------------------------------------- */
/*
std::vector<unsigned int> Moves::detect_molecule(
    std::vector<std::vector<unsigned int> > neigh_list,
    const std::vector<Particle> &particles, std::vector<Particle> molecule,
     bool &detected)
{
    unsigned int i, count;
    std::vector<unsigned int> elm_idx;
    
    count = 0;
    while (count < neigh_list.size() && !detected) {
        elm_idx.clear();
        i = rng->next_int(neigh_list.size());     // pick initial particle
        check_neigh_recu(i, particles, molecule, 0, elm_idx, neigh_list);
        if (elm_idx.size() == molecule.size()) {
            detected = true;
        }
        count ++;
    }
    if (!detected) {
        elm_idx.clear();
    }
    return elm_idx;
}
*/
