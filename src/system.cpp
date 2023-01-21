/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-23)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  This where the system class is implemented, with all its functionality. The
  system may contain several boxes. The moves are associated with the system
  class, as a move may touch several boxes. Most other properties (boundary,
  forcefield, constraints and so on) are associated with a box. System class
  includes all keyword initialization methods.
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>
#include <functional>
#include <map>
#include <cmath>

#include "system.h"
#include "box.h"
#include "tqdm.h"
#include "dump.h"
#include "thermo.h"
#include "particle.h"
#include "distance_manager.h"

#include "rng/rng.h"
#include "rng/simple.h"
#include "rng/mersennetwister.h"
#include "rng/fastmersennetwister.h"

#include "boundary/boundary.h"
#include "boundary/open.h"
#include "boundary/periodic.h"

#include "forcefield/forcefield.h"
#include "forcefield/idealgas.h"
#include "forcefield/lennardjones.h"
#include "forcefield/vashishta.h"

#include "sampler/sampler.h"
#include "sampler/metropolis.h"
#include "sampler/umbrella.h"

#include "moves/moves.h"
#include "moves/trans.h"
#include "moves/transmh.h"
#include "moves/avbmc.h"
#include "moves/avbmcin.h"
#include "moves/avbmcout.h"
#include "moves/avbmcmol.h"
#include "moves/avbmcmolin.h"
#include "moves/avbmcmolout.h"
#include "moves/avbmcswapright.h"
#include "moves/avbmcmolswapright.h"

#include "constraint/constraint.h"
#include "constraint/maxdistance.h"
#include "constraint/mindistance.h"
#include "constraint/maxneigh.h"
#include "constraint/minneigh.h"


/* ----------------------------------------------------------------------------
   System constructor, taking the working directory 'working_dir_in' as
   argument. The argument 'memory_intensity' specifies how memeory intensive
   the simulation should be, with a memory-cpu-time tradeoff. The options are:

       1: Storing necessary neighbor lists only
       2: Storing distances and relative coordinates between particles
       3: Storing distances, relative coordinates and energy contributions
          of each particle
---------------------------------------------------------------------------- */

System::System(const std::string &working_dir_in, bool initialize_in) 
    : working_dir(working_dir_in)
{
    nbox = nmove = step = 0;
    time = temp = chempot = poteng = beta = 0.;
    ndim = 3;

    //print_logo();
    logo_printed = false;
    rng_allocated_externally = false;
    sampler_allocated_externally = false;

    // set default rng and sampler
    rng = new MersenneTwister;
    sampler = new Metropolis(this);

    // initialize box 
    if (initialize_in) {
        Box *box = new Box(this);
        add_box(box);
        box->box_allocated_in_system = true;
    }
}


/* ----------------------------------------------------------------------------
   Copy constructor
---------------------------------------------------------------------------- */

System::System(const System &other) 
    : working_dir(other.working_dir), boxes(other.boxes), moves(other.moves),
      moves_prob(other.moves_prob)
{
    time = other.time;
    temp = other.temp;
    //logo_printed = other.logo_printed;
    chempot = other.chempot;
    poteng = other.poteng;
    nbox = other.nbox;
    nmove = other.nmove;
    step = other.step;
    ndim = other.ndim;
    rng_allocated_externally = other.rng_allocated_externally;
    sampler_allocated_externally = other.sampler_allocated_externally;
}


/* ----------------------------------------------------------------------------
   Make and move to working directory
---------------------------------------------------------------------------- */

void System::set_working_directory(const std::string &working_directory)
{
    //mkdir(working_dir);
    //chdir(working_dir);
}


/* ----------------------------------------------------------------------------
   Set number of dimensions
---------------------------------------------------------------------------- */

void System::set_dim(unsigned int dim)
{
    for (Box *box : boxes) {
        if (box->npar > 0) {
            std::cout << "Cannot change dimensionality when particles are "
                      << "added! Aborting." << std::endl;
            exit(0);
        }
    }
    ndim = dim;
}


/* ----------------------------------------------------------------------------
   Set box temperature used in NVT and uPT Monte Carlo simulations. The 
   Boltzmann constant can also be set to match the units. Default kB is
   1., which works for many scaled units (for instance Lennard-Jones units)
---------------------------------------------------------------------------- */

void System::set_temp(const double temp_in, const double kB_in)
{
    temp = temp_in;
    kB = kB_in;
    beta = 1. / (kB * temp);
}


/* ----------------------------------------------------------------------------
   Set chemical potential of system used in grand canonical ensemble. Instead
   of setting the chemical potential directly, it can be indirectly set by
   the saturation, which is ensemble invariant. This should be done after
   temperature and the Boltzmann constant are set.
---------------------------------------------------------------------------- */

void System::set_chempot(const double chempot_in, const double sat_in)
{
    if (sat_in > 0) {
        chempot = - std::log(sat_in) / beta;
    }
    else {
        chempot = chempot_in;
    }
}


/* ----------------------------------------------------------------------------
   Set seed to be used by RNG
---------------------------------------------------------------------------- */

void System::set_seed(unsigned int seed_)
{
    rng->set_seed(seed_);
}


/* ----------------------------------------------------------------------------
   Set mass of chemical symbol. Masses of all chemical symbols have to be given
   if running molecular dynamics simulations, as the software does not look up
   the masses in a table.
---------------------------------------------------------------------------- */
/*
void System::set_mass(const std::string label, const double mass)
{
    mass_labels.push_back(label);
    masses.push_back(mass);
}
*/


/* ----------------------------------------------------------------------------
   Overwrite default sampler, Metropolis
---------------------------------------------------------------------------- */

void System::set_sampler(class Sampler *sampler_in)
{
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    sampler = sampler_in;
    sampler_allocated_externally = true;
}


void System::set_sampler(const std::string &sampler_in,
    std::function<double(int)> f, int ntabulated)
{
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    if (sampler_in == "metropolis") {
        sampler = new Metropolis(this);
        sampler_allocated_externally = false;
    }
    else if (sampler_in == "umbrella") {
        sampler = new Umbrella(this, f, ntabulated);
        sampler_allocated_externally = false;
    }
    else {
        std::cout << "Sampler '" << sampler_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


void System::set_sampler(const std::string &sampler_in,
    std::valarray<double> tabulated)
{
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    if (sampler_in == "metropolis") {
        sampler = new Metropolis(this);
        sampler_allocated_externally = false;
    }
    else if (sampler_in == "umbrella") {
        sampler = new Umbrella(this, tabulated);
        sampler_allocated_externally = false;
    }
    else {
        std::cout << "Sampler '" << sampler_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Overwrite default random number generator, Mersenne Twister
---------------------------------------------------------------------------- */

void System::set_rng(class RandomNumberGenerator* rng_in)
{
    if (!rng_allocated_externally) {
        delete rng;
    }
    rng = rng_in;
    rng_allocated_externally = true;
}


void System::set_rng(const std::string &rng_in)
{
    if (!rng_allocated_externally) {
        delete rng;
    }
    if (rng_in == "simple") {
        rng = new Simple;
        rng_allocated_externally = false;
    }
    else if (rng_in == "mersennetwister") {
        rng = new MersenneTwister;
        rng_allocated_externally = false;
    }
    else if (rng_in == "mt19937") {
        rng = new MersenneTwister;
        rng_allocated_externally = false;
    }
    else if (rng_in == "fastmersennetwister") {
        rng = new FastMersenneTwister;
        rng_allocated_externally = false;
    }
    else {
        std::cout << "Random number generator '" << rng_in 
                  << "' is not implemented! Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Set forcefield. This method is forward to the box and can only be done if a
   box was detected.
---------------------------------------------------------------------------- */

void System::set_forcefield(class ForceField* forcefield_in,
    int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! "
                  << "Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            if (box->forcefield_allocated_in_system) {
                delete box->forcefield;
            }
            box->set_forcefield(forcefield_in);
            box->forcefield_allocated_in_system = false;
        }
    }
    else {
        if (boxes[box_id]->forcefield_allocated_in_system) {
            delete boxes[box_id]->forcefield;
        }
        boxes[box_id]->set_forcefield(forcefield_in);
        boxes[box_id]->forcefield_allocated_in_system = false;
    }
}


void System::set_forcefield(const std::string &forcefield_in,
    const std::vector<std::string> &labels, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! "
                  << "Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (forcefield_in == "idealgas") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                if (box->forcefield_allocated_in_system) {
                    delete box->forcefield;
                }
                box->forcefield = new IdealGas(box, labels);
                box->forcefield_allocated_in_system = true;
                box->initialized = true;
            }
        }
        else {
            if (boxes[box_id]->forcefield_allocated_in_system) {
                delete boxes[box_id]->forcefield;
            }
            boxes[box_id]->forcefield = new IdealGas(boxes[box_id], labels);
            boxes[box_id]->forcefield_allocated_in_system = true;
            boxes[box_id]->initialized = true;
        }
    }
    else {
        std::cout << "Force-field '" << forcefield_in << "' is not implemented"
                  << "or does not have the given signature! Aborting."
                  << std::endl;
        exit(0);
    }
}


void System::set_forcefield(const std::string &forcefield_in,
    const std::string &paramfile, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (forcefield_in == "lennardjones") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                if (box->forcefield_allocated_in_system) {
                    delete box->forcefield;
                }
                box->forcefield = new LennardJones(box, paramfile);
                box->forcefield_allocated_in_system = true;
                box->initialized = true;
            }
        }
        else {
            if (boxes[box_id]->forcefield_allocated_in_system) {
                delete boxes[box_id]->forcefield;
            }
            boxes[box_id]->forcefield = new LennardJones(boxes[box_id], paramfile);
            boxes[box_id]->forcefield_allocated_in_system = true;
            boxes[box_id]->initialized = true;
        }
    }
    else if (forcefield_in == "vashishta") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                if (box->forcefield_allocated_in_system) {
                    delete box->forcefield;
                }
                box->forcefield = new Vashishta(box, paramfile);
                box->forcefield_allocated_in_system = true;
                box->initialized = true;
            }
        }
        else {
            if (boxes[box_id]->forcefield_allocated_in_system) {
                delete boxes[box_id]->forcefield;
            }
            boxes[box_id]->forcefield = new Vashishta(boxes[box_id], paramfile);
            boxes[box_id]->forcefield_allocated_in_system = true;
            boxes[box_id]->initialized = true;
        }
    }
    else {
        std::cout << "Force-field '" << forcefield_in << "' is not implemented or"
                  << "does not have the given signature! Aborting." << std::endl;
        exit(0);
    }
}
   

/* ----------------------------------------------------------------------------
   Set boundary. This method is forward to the box and can only be done if a
   box was detected.
---------------------------------------------------------------------------- */

void System::set_boundary(class Boundary *boundary_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set boundary" << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            if (!box->boundary_allocated_externally) {
                delete box->boundary;
                box->boundary_allocated_externally = true;
            }
            if (box->boundary_allocated_in_system) {
                delete box->boundary;
            }
            box->set_boundary(boundary_in);
            box->boundary_allocated_in_system = false;
        }
    }
    else {
        boxes[box_id]->set_boundary(boundary_in);
        boxes[box_id]->boundary_allocated_in_system = false;
    }
}


void System::set_boundary(const std::string &boundary_in,
    std::valarray<double> length, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set boundary" << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (boundary_in == "open") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                if (!box->boundary_allocated_externally) {
                    delete box->boundary;
                    box->boundary_allocated_externally = true;
                }
                if (box->boundary_allocated_in_system) {
                    delete box->boundary;
                }
                box->boundary = new Open(box);
                box->boundary_allocated_in_system = true;
            }
        }
        else {
            if (!boxes[box_id]->boundary_allocated_externally) {
                delete boxes[box_id]->boundary;
                boxes[box_id]->boundary_allocated_externally = true;
            }
            if (boxes[box_id]->boundary_allocated_in_system) {
                delete boxes[box_id]->boundary;
            }
            boxes[box_id]->boundary = new Open(boxes[box_id]);
            boxes[box_id]->boundary_allocated_in_system = true;
        }
    }
    else if (boundary_in == "periodic") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                if (!box->boundary_allocated_externally) {
                    delete box->boundary;
                    box->boundary_allocated_externally = true;
                }
                if (box->boundary_allocated_in_system) {
                    delete box->boundary;
                }
                box->boundary = new Periodic(box, length);
                box->boundary_allocated_in_system = true;
            }
        }
        else {
            if (!boxes[box_id]->boundary_allocated_externally) {
                delete boxes[box_id]->boundary;
                boxes[box_id]->boundary_allocated_externally = true;
            }
            if (boxes[box_id]->boundary_allocated_in_system) {
                delete boxes[box_id]->boundary;
            }
            boxes[box_id]->boundary = new Periodic(boxes[box_id], length);
            boxes[box_id]->boundary_allocated_in_system = true;
        }
    }
    else {
        std::cout << "Boundary '" << boundary_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Add constraint. This method is forward to the box and can only be done if a
   box was detected.
---------------------------------------------------------------------------- */

void System::add_constraint(Constraint* constraint_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add constraint" << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        if (nbox > 1) {
            std::cout << "Warning: More than one box was detected. Setting "
                      << "constraint for all boxes" << std::endl;
        }
        for (Box *box : boxes) {
            box->add_constraint(constraint_in);
        }
    }
    else {
        boxes[box_id]->add_constraint(constraint_in);
    }
}


void System::add_constraint(const std::string &constraint_in,
    const std::string &element1, const std::string &element2, double distance,
    int nneigh, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add constraint" << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }
    else if (box_id < 0 && nbox > 1) {
        std::cout << "Warning: More than one box was detected. Setting "
                  << "constraint for all boxes" << std::endl;
    }

    if (constraint_in == "maxdistance") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Constraint *constraint = new MaxDistance(box, element1, element2, distance);
                box->add_constraint(constraint);
                box->constraint_allocated_in_system[box->nconstraint-1] = true;
            }
        }
        else {
            Constraint *constraint = new MaxDistance(boxes[box_id], element1, element2, distance);
            boxes[box_id]->add_constraint(constraint);
            boxes[box_id]->constraint_allocated_in_system[boxes[box_id]->nconstraint-1] = true;
        }
    }
    else if (constraint_in == "mindistance") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Constraint *constraint = new MinDistance(box, element1, element2, distance);
                box->add_constraint(constraint);
                box->constraint_allocated_in_system[box->nconstraint-1] = true;
            }
        }
        else {
            Constraint *constraint = new MinDistance(boxes[box_id], element1, element2, distance);
            boxes[box_id]->add_constraint(constraint);
            boxes[box_id]->constraint_allocated_in_system[boxes[box_id]->nconstraint-1] = true;
        }
    }
    else if (constraint_in == "maxneigh") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Constraint *constraint = new MaxNeigh(box, element1, element2, distance, nneigh);
                box->add_constraint(constraint);
                box->constraint_allocated_in_system[box->nconstraint-1] = true;
            }
        }
        else {
            Constraint *constraint = new MaxNeigh(boxes[box_id], element1, element2, distance, nneigh);
            boxes[box_id]->add_constraint(constraint);
            boxes[box_id]->constraint_allocated_in_system[boxes[box_id]->nconstraint-1] = true;
        }
    }
    else if (constraint_in == "minneigh") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Constraint *constraint = new MinNeigh(box, element1, element2, distance, nneigh);
                box->add_constraint(constraint);
                box->constraint_allocated_in_system[box->nconstraint-1] = true;
            }
        }
        else {
            Constraint *constraint = new MinNeigh(boxes[box_id], element1, element2, distance, nneigh);
            boxes[box_id]->add_constraint(constraint);
            boxes[box_id]->constraint_allocated_in_system[boxes[box_id]->nconstraint-1] = true;
        }
    }
    else if (constraint_in == "stillinger") {
        std::cout << "Stillinger constraint cannot be initialized this way! \n"
                  << "Use object method: stillinger = avbmc.Stillinger(box) \n"
                  << "stillinger.set_criterion(element1, element2, distance) " << std::endl;
        exit(0);
    }
    else {
        std::cout << "Constraint '" << constraint_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Remove constraint of box 'box_id' by index
---------------------------------------------------------------------------- */

void System::rm_constraint(unsigned int idx, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    boxes[box_id]->rm_constraint(idx);
}


/* ----------------------------------------------------------------------------
   Take a snapshot of a given box.
---------------------------------------------------------------------------- */

void System::snapshot(const std::string &filename, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot take snapshot! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    boxes[box_id]->snapshot(filename);
}


/* ----------------------------------------------------------------------------
   Set dump. This method is forward to the box and can only be done if
   a box was detected.
---------------------------------------------------------------------------- */

void System::set_dump(int freq, const std::string &filename,
    const std::vector<std::string> &outputs, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set dump! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    boxes[box_id]->set_dump(freq, filename, outputs);
}


/* ----------------------------------------------------------------------------
   Set thermo. This method is forward to the box and can only be done if a box
   was detected.
---------------------------------------------------------------------------- */

void System::set_thermo(int freq, const std::string &filename,
    const std::vector<std::string> &outputs, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set thermo! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    boxes[box_id]->set_thermo(freq, filename, outputs);
}


/* ----------------------------------------------------------------------------
   Add move type and the corresponding probability. The probabilities have to
   add up to 1. Forcefield has to initialized first, to link AVBMC atoms to
   types.
---------------------------------------------------------------------------- */

void System::add_move(Moves* move, double prob)
{
    /*
    if (!initialized) {
        std::cout << "Forcefield needs to be initialized before adding moves!"
                  << std::endl;
        exit(0);
    }
    */
    nmove ++;
    moves.push_back(move);
    moves_prob.push_back(prob);
    moves_allocated_in_system.push_back(false);
}


void System::add_move(const std::string &move_in, double prob, double dx,
    double Ddt, const std::string &element, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (move_in == "trans") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new Trans(this, box, dx, element);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new Trans(this, boxes[box_id], dx, element);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else if (move_in == "transmh") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new TransMH(this, box, dx, Ddt, element);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new TransMH(this, boxes[box_id], dx, Ddt, element);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else {
        std::cout << "Move '" << move_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


void System::add_move(const std::string &move_in, double prob, 
    const std::string &particle_in, double r_below, double r_above,
    bool energy_bias, int box_id, int box_id2)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (move_in == "avbmc") {
        // simply add both avbmcin and avbmcout with equal probability
        add_move("avbmcin", prob / 2., particle_in, r_below, r_above, energy_bias, box_id);
        add_move("avbmcout", prob / 2., particle_in, r_below, r_above, energy_bias, box_id);
    }
    else if (move_in == "avbmcin") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new AVBMCIn(this, box, particle_in, r_below, r_above, energy_bias);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new AVBMCIn(this, boxes[box_id], particle_in, r_below, r_above, energy_bias);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else if (move_in == "avbmcout") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new AVBMCOut(this, box, particle_in, r_above, energy_bias);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new AVBMCOut(this, boxes[box_id], particle_in, r_above, energy_bias);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else if (move_in == "avbmcswapright") {
        if (box_id < 0 || box_id2 < 0) {
            std::cout << "Both box_id1 and box_id2 have to be defined in order to do "
                      << "inter-box swap moves! Aborting." << std::endl;
            exit(0);
        }
        else if (box_id2 >= nbox)
        {
            std::cout << "Box-ID 2 is out of range! " << nbox << " boxes found. "
                      << "Aborting." << std::endl;
            exit(0);
        }

        Moves *move = new AVBMCSwapRight(this, boxes[box_id], boxes[box_id2], particle_in, r_below, r_above, energy_bias);
        add_move(move, prob);
        moves_allocated_in_system[nmove-1] = true;
    }
    else if (move_in == "avbmcswap") {
        add_move("avbmcswapright", prob / 2., particle_in, r_below, r_above, energy_bias, box_id, box_id2);
        add_move("avbmcswapright", prob / 2., particle_in, r_below, r_above, energy_bias, box_id2, box_id);
    }
    else {
        std::cout << "Move '" << move_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


void System::add_move(const std::string &move_in, double prob,
    std::vector<Particle> molecule_in, double r_below,
    double r_above, double r_inner, bool energy_bias, bool target_mol,
    int box_id, int box_id2)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (move_in == "avbmcmol") {
        // simply add both avbmcmolin and avbmcmolout with equal probability
        add_move("avbmcmolin", prob / 2., molecule_in, r_below, r_above, r_inner, energy_bias, target_mol, box_id);
        add_move("avbmcmolout", prob / 2., molecule_in, r_below, r_above, r_inner, energy_bias, target_mol, box_id);
    }
    else if (move_in == "avbmcmolin") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new AVBMCMolIn(this, box, molecule_in, r_below, r_above, r_inner, energy_bias, target_mol);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new AVBMCMolIn(this, boxes[box_id], molecule_in, r_below, r_above, r_inner, energy_bias, target_mol);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else if (move_in == "avbmcmolout") {
        if (box_id < 0) {
            for (Box *box : boxes) {
                Moves *move = new AVBMCMolOut(this, box, molecule_in, r_above, r_inner, energy_bias, target_mol);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else {
            Moves *move = new AVBMCMolOut(this, boxes[box_id], molecule_in, r_above, r_inner, energy_bias, target_mol);
            add_move(move, prob);
            moves_allocated_in_system[nmove-1] = true;
        }
    }
    else if (move_in == "avbmcmolswapright") {
        if (box_id < 0 || box_id2 < 0) {
            std::cout << "Both box_id1 and box_id2 have to be defined in order to do "
                      << "inter-box swap moves! Aborting." << std::endl;
            exit(0);
        }
        else if (box_id2 >= nbox)
        {
            std::cout << "Box-ID 2 is out of range! " << nbox << " boxes found. "
                      << "Aborting." << std::endl;
            exit(0);
        }

        Moves *move = new AVBMCMolSwapRight(this, boxes[box_id], boxes[box_id2], molecule_in, r_below, r_above, r_inner, energy_bias, target_mol);
        add_move(move, prob);
        moves_allocated_in_system[nmove-1] = true;
    }
    else if (move_in == "avbmcmolswap") {
        add_move("avbmcmolswapright", prob / 2., molecule_in, r_below, r_above, r_inner, energy_bias, target_mol, box_id, box_id2);
        add_move("avbmcmolswapright", prob / 2., molecule_in, r_below, r_above, r_inner, energy_bias, target_mol, box_id2, box_id);
    }
    else {
        std::cout << "Move '" << move_in << "' is not implemented! "
                  << "Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Remove move by index 'idx'
---------------------------------------------------------------------------- */

void System::rm_move(unsigned int idx)
{
    if (moves_allocated_in_system[idx]) {
        delete moves[idx];
    }
    moves.erase(moves.begin() + idx);
    moves_prob.erase(moves_prob.begin() + idx);
    moves_allocated_in_system.erase(moves_allocated_in_system.begin() + idx);
    nmove--;
}


/* ----------------------------------------------------------------------------
   Add box 'box_in' to system
---------------------------------------------------------------------------- */

void System::add_box(Box* box_in)
{
    box_in->box_id = nbox;
    boxes.push_back(box_in);
    nbox ++;
}


void System::add_box()
{
    Box *box = new Box(this);
    add_box(box);
    box->box_allocated_in_system = true;
}


/* ----------------------------------------------------------------------------
   Remove box by index 'box_id'. All later box indices will be shifted
---------------------------------------------------------------------------- */

void System::rm_box(unsigned int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (boxes[box_id]->forcefield_allocated_in_system) {
        delete boxes[box_id]->forcefield;
    }
    if (boxes[box_id]->boundary_allocated_in_system) {
        delete boxes[box_id]->boundary;
    }
    for (unsigned int i=boxes[box_id]->nconstraint; i--;) {
        rm_constraint(i, box_id=box_id);
    }
    /*
    for (unsigned int i=0; i<boxes[box_id]->nconstraint; i++) {
        if (boxes[box_id]->constraint_allocated_in_system[i]) {
            delete boxes[box_id]->constraints[i];
        }
    }
    */
    if (boxes[box_id]->box_allocated_in_system) {
        delete boxes[box_id];
    }
    boxes.erase(boxes.begin() + box_id);
    nbox--;
}


/* ----------------------------------------------------------------------------
   Add particles to the system. 
---------------------------------------------------------------------------- */

void System::add_particle(Particle particle_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            box->add_particle(particle_in);
        }
    }
    else {
        boxes[box_id]->add_particle(particle_in);
    }
}


void System::add_particle(const std::string &label_in,
    std::valarray<double> pos, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            box->add_particle(label_in, pos);
        }
    }
    else {
        boxes[box_id]->add_particle(label_in, pos);
    }
}


void System::add_particles(std::vector<Particle> particles_in,
    int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            box->add_particles(particles_in);
        }
    }
    else {
        boxes[box_id]->add_particles(particles_in);
    }
}


void System::add_particles(const std::string &label_in,
    std::vector<std::valarray<double> > pos, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            box->add_particles(label_in, pos);
        }
    }
    else {
        boxes[box_id]->add_particles(label_in, pos);
    }
}


/* ----------------------------------------------------------------------------
   Read particles from an xyz-file
---------------------------------------------------------------------------- */

void System::read_particles(const std::string &filename, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    if (box_id < 0) {
        for (Box *box : boxes) {
            box->read_particles(filename);
        }
    }

    boxes[box_id]->read_particles(filename);
}


/* ----------------------------------------------------------------------------
   Remove particle from box 'box_id' by index 'idx'
---------------------------------------------------------------------------- */

void System::rm_particle(unsigned int idx, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot remove particles! Aborting."
                  << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }

    boxes[box_id]->rm_particle(idx);
}


void System::clear_particles(int box_id)
{
    if (box_id < 0) {
        for (Box *box : boxes) {
            box->clear_particles();
        }
    }
    else {
        boxes[box_id]->clear_particles();
    }
}


/* ----------------------------------------------------------------------------
   Write histogram of system sizes out to file
---------------------------------------------------------------------------- */

void System::write_size_histogram(const std::string &filename,
    int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }
    boxes[box_id]->write_size_histogram(filename);
}


std::vector<unsigned int> System::get_size_histogram(int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found! Aborting." << std::endl;
        exit(0);
    }
    else if (box_id >= nbox)
    {
        std::cout << "Box-ID is out of range! " << nbox << " boxes found. "
                  << "Aborting." << std::endl;
        exit(0);
    }
    return boxes[box_id]->size_histogram;
}

/* ----------------------------------------------------------------------------
   Returns the last iteration
---------------------------------------------------------------------------- */

unsigned int System::get_maxiter(unsigned int nsteps)
{
    unsigned int maxiter;

    maxiter = 0;
    if (step == 0) {
        maxiter += 1;
    }
    return maxiter + step + nsteps;
}


/* ----------------------------------------------------------------------------
   Print logo header
---------------------------------------------------------------------------- */

void System::print_logo()
{
    std::cout << std::endl;
    std::cout << " █████╗ ██╗   ██╗██████╗ ███╗   ███╗ ██████╗" << std::endl;
    std::cout << "██╔══██╗██║   ██║██╔══██╗████╗ ████║██╔════╝" << std::endl;
    std::cout << "███████║██║   ██║██████╔╝██╔████╔██║██║     " << std::endl;
    std::cout << "██╔══██║╚██╗ ██╔╝██╔══██╗██║╚██╔╝██║██║     " << std::endl;
    std::cout << "██║  ██║ ╚████╔╝ ██████╔╝██║ ╚═╝ ██║╚██████╗" << std::endl;
    std::cout << "╚═╝  ╚═╝  ╚═══╝  ╚═════╝ ╚═╝     ╚═╝ ╚═════╝" << std::endl;
    std::cout << std::endl;
    logo_printed = true;
}


/* ----------------------------------------------------------------------------
   Print information about simulation
---------------------------------------------------------------------------- */

void System::print_info()
{
    std::time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Started computation on " << std::ctime(&start_time) << std::endl;

    //std::cout << std::fixed;
    //std::cout << std::boolalpha;
    //std::cout << std::setprecision(6);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "            System Information " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "Number of dimensions:     " << ndim << std::endl;
    //std::cout << "Forcefield:               " << forcefield->label 
    //          << " " << forcefield->paramfile << std::endl;
    std::cout << "Random number generator:  " << rng->label << std::endl;
    std::cout << std::endl;
    //std::cout << "Number of particle types: " << forcefield->ntype << std::endl;
    //std::cout << "Unique particle types:";
    //for(int i=0; i < forcefield->ntype; i++){
    //    std::cout << " " << forcefield->unique_labels[i];
    //}
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Number of boxes: " << nbox << std::endl;
    for(unsigned int i=0; i < nbox; i++){
        std::cout << "  Box " << i+1 << ":" << std::endl;
        std::cout << "    Number of atoms:   " << boxes[i]->npar << std::endl;
        std::cout << "    Boundary:          " << boxes[i]->boundary->label << std::endl;
    }
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
   Print Monte Carlo information
------------------------------------------------------------------------------- */

void System::print_mc_info()
{
    std::cout << std::endl;
    std::cout << "           Monte Carlo Information " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "Temperature:              " << temp << std::endl;
    std::cout << "Chemical potential:       " << chempot << std::endl;
    std::cout << "Sampling method:          " << sampler->label << std::endl;
    std::cout << std::endl;
    std::cout << "Number of move types:     " << nmove << std::endl;
    for(unsigned int i=0; i < nmove; i++){
        std::cout << "  Move type " << i+1 << ": ";
        std::cout << moves[i]->repr();
        std::cout << "    Probability: " << moves_prob[i] << std::endl;
    }
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
   Run molecular dynamics simulation
---------------------------------------------------------------------------- */
/*
void Box::run_md(unsigned int nsteps)
{
    init_simulation();
    check_masses();
    print_info();
    thermo->print_header();

    int maxiter = get_maxiter(nsteps);

    // run molecular dynamics simulation
    auto t1 = chrono::high_resolution_clock::now();
    while(step<maxiter){
        time = step * integrator->dt;
        dump->print_frame();
        thermo->print_line();
        energy = integrator->next_step();
        step ++;
    }
    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = chrono::duration<double>(t2 - t1).count();
    cout << "Elapsed time: " << duration_seconds << "s" << endl;
}
*/

/* ----------------------------------------------------------------------------
   Run Monte Carlo simulation with 'nsteps' cycles and 'nmoves' moves per
   cycle. A progress bar is printed using tqdm. 
---------------------------------------------------------------------------- */

void System::initialize_mc_run()
{
    for (Box* box : boxes) {
        box->size_histogram.resize(box->npar + 1);
        box->size_histogram[box->npar] ++;
        if (box->store_energy) {
            box->forcefield->initialize();
        }
    }

    for (Moves* move : moves) {
        move->ndrawn = 0;
        move->naccept = 0;
    }

    double sum_prob = std::accumulate(moves_prob.begin(), moves_prob.end(), 0.);
    assert ((sum_prob - 1.0) < 0.01);
}


/* ----------------------------------------------------------------------------
   Run nsteps cycles
---------------------------------------------------------------------------- */

void System::run_mc(const unsigned int nsteps, const unsigned int nmoves)
{
    initialize_mc_run();

    unsigned int init_step = step;
    tqdm bar;
    unsigned int maxiter = get_maxiter(nsteps);
    while (step < maxiter) {
        bar.progress(step-init_step, nsteps);
        run_mc_cycle(nmoves);
    }
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
   Run a cycle
---------------------------------------------------------------------------- */

void System::run_mc_cycle(const unsigned int nmoves)
{
    for (Box* box : boxes) {
        box->dump->print_frame(step);
        box->thermo->print_line(step);
    }
    sampler->sample(nmoves);
    step ++;
}


/* ----------------------------------------------------------------------------
   Print beautiful table for moves statistics
---------------------------------------------------------------------------- */

std::string cell_padding(const std::string &text, const std::string &side, std::size_t column_width)
{
    std::size_t whitespace, whitespace_left, whitespace_right;

    if (text.size() > column_width) return text;
    
    whitespace = column_width - text.size();
    if (side == "center" || side == "c") {
        whitespace_left = static_cast<int>(std::floor(whitespace/2.));
        whitespace_right = static_cast<int>(std::ceil(whitespace/2.));
    }
    else if (side == "left" || side == "l") {
        whitespace_left = 0;
        whitespace_right = whitespace;
    }
    else if (side == "right" || side == "r") {
        whitespace_left = whitespace;
        whitespace_right = 0;
    }
    else {
        exit(0);
    }

    return std::string(whitespace_left, ' ') + text + std::string(whitespace_right, ' ');
}


std::string create_table(std::vector<std::string> keywords,
    std::vector<std::vector<std::string> > values, const std::string &style = "basic", bool print = true) 
{
    std::string table, pos_col;

    pos_col = "l" + std::string(keywords.size(), 'r');

    if (style == "latex") {
        table = "\\begin{center}\n";
        std::cout << table << std::endl;
        table += "  \\begin{tabular}{" + pos_col + "}\\hline\n";
        std::cout << table << std::endl;
        table += "    ";
        std::cout << table << std::endl;
        // keyword
        for (std::string keyword : keywords) {
            table += keyword + " & ";
        }
        std::cout << table << std::endl;
        table += " \\\\ \n";
        std::cout << table << std::endl;
        table += "    \\hline\\hline\n";
        std::cout << table << std::endl;
        // values
        for (std::vector<std::string> row : values) {
            table += "    ";
            for (std::string value : row) {
                std::cout << value << std::endl;
                table += value + " & ";
            }
            table += " \\\\ \n";
        }
        table += "    \\hline \n";
        table += "  \\end{tabular} \n";
        table += "\\end{center} \n";
    }
    else if (style == "basic") {
        std::string head_row, row, rows, total;

        std::valarray<std::size_t> max_size_col(keywords.size());
        std::valarray<double> sums(0., keywords.size());
        std::valarray<std::string> sums_str(keywords.size());

        // determine max size of column
        for (std::size_t i=0; i<keywords.size(); i++) {
            max_size_col[i] = keywords[i].size();
        }
        for (std::vector<std::string> value : values) { 
            for (std::size_t i=0; i<value.size(); i++) {
                std::size_t len_cell = value[i].size();
                try {
                    sums[i] += std::stod(value[i]);
                    sums_str[i] = std::to_string(sums[i]);
                }
                catch (...) {
                    sums_str[i] = "";
                }
                std::size_t len_tot = sums_str[i].size();
                if (len_cell > max_size_col[i]) max_size_col[i] = len_cell;
                if (len_tot > max_size_col[i]) max_size_col[i] = len_tot;
            }
        }

        // get header
        for (std::size_t i=0; i<keywords.size(); i++) {
            std::string cell = keywords[i];
            head_row += "| " + cell_padding(cell, "c", max_size_col[i]) + " ";
        }
        head_row += "|";

        // get rows
        for (std::vector<std::string> value : values) {
            row = "";
            for (std::size_t i=0; i<value.size(); i++) {
                row += "| " + cell_padding(value[i], std::to_string(pos_col[i]), max_size_col[i]) + " ";
            }
            row += "|";
            rows += row + "\n";
        }

        // total
        for (std::size_t i=0; i<sums_str.size(); i++) {
            total += "| " + cell_padding(sums_str[i], std::to_string(pos_col[i]), max_size_col[i]) + " ";
        }
        total += "|\n";

        // merging everything together
        table += std::string(head_row.size(), '-') + "\n"; 
        table += head_row + "\n";
        table += std::string(head_row.size(), '-') + "\n";
        table += rows;
        table += std::string(head_row.size(), '-') + "\n"; 
        table += total;
        table += std::string(head_row.size(), '-') + "\n";
    }

    if (print) {
        std::cout << table << std::endl;
    }

    return table;
}


std::string get_column_mapping(Moves *move, const std::string &keyword)
{
    std::string rejout, rejtarg, rejtargout, rejtargin;
    std::vector<std::string> mylist{"AVBMCMolOut", "AVBMCOut", "AVBMCSwapRight", "AVBMCMolSwapRight"};

    if (std::find(std::begin(mylist), std::end(mylist), move->label) != std::end(mylist)) {
        rejout = std::to_string(move->nrejectout);
        rejtarg = std::to_string(move->nrejecttarget);
    }
    else if (move->label == "AVBMCMolIn") {
        rejout = "-";
        rejtarg = std::to_string(move->nrejecttarget);
    }
    else {
        rejout = rejtarg = "-";
    }

    if (move->label == "AVBMCMolSwapRight") {
        rejtargout = std::to_string(move->nrejecttargetout);
        rejtargin = std::to_string(move->nrejecttargetin);
    }
    else {
        rejtargout = rejtargin = "-";
    }

    std::map<std::string, std::string> column_map { 
        {"move", move->label},
        {"ndrawn", std::to_string(move->ndrawn)},
        {"naccept", std::to_string(move->naccept)},
        {"nreject", std::to_string(move->ndrawn-move->naccept)},
        {"accratio", std::to_string(static_cast<double>(move->naccept)/move->ndrawn)},
        {"cputime", std::to_string(move->cum_time)},
        {"rejout", rejout},
        {"rejtarg", rejtarg},
        {"rejtargout", rejtargout},
        {"rejtargin", rejtargin}
    };

    return column_map[keyword];
}


std::string System::print_statistics(std::vector<std::string> cols, bool print, const std::string &style)
{
    std::string head_row, row, rows, total, table, pad;

    std::map<std::string, std::string> head_labels {
        {"move", "Move"},
        {"ndrawn", "#drawn"},
        {"naccept", "#accept"},
        {"nreject", "#reject"},
        {"accratio", "acc. ratio"},
        {"cputime", "CPU-time (s)"},
        {"rejout", "#reject out"},
        {"rejtarg", "#reject target"},
        {"rejtargout", "#reject target out"},
        {"rejtargin", "#reject target in"}
    };

    std::vector<std::size_t> max_size_col(cols.size());

    // determine max size of column
    for (std::size_t i=0; i<cols.size(); i++) {
        max_size_col[i] = head_labels[cols[i]].size();
        for (Moves *move : moves) {
            std::size_t len_cell = get_column_mapping(move, cols[i]).size();
            if (len_cell > max_size_col[i]) max_size_col[i] = len_cell;
        }
    }

    // get header
    for (std::size_t i=0; i<cols.size(); i++) {
        pad = "c";
        std::string cell = head_labels[cols[i]];
        head_row += "| " + cell_padding(cell, pad, max_size_col[i]) + " ";
    }
    head_row += "|";

    // get rows
    std::valarray<double> totals(0., cols.size());
    for (Moves *move : moves) {
        row = "";
        for (std::size_t i=0; i<cols.size(); i++) {
            std::string cell = get_column_mapping(move, cols[i]);
            if (cols[i] == "move") {
                pad = "l";
            }
            else {
                pad = "r";
                if (cell != "-") {
                    totals[i] += std::stod(cell);
                }
            }
            row += "| " + cell_padding(cell, pad, max_size_col[i]) + " ";
        }
        row += "|";
        rows += row + "\n";
    }

    // get total
    for (std::size_t i=0; i<cols.size(); i++) {
        if (cols[i] == "move") {
            total += "| " + cell_padding("Total:", "l", max_size_col[i]) + " ";
        }
        else if (cols[i] == "cputime") {
            total += "| " + cell_padding(std::to_string(totals[i]), "r", max_size_col[i]) + " ";
        }
        else if (cols[i] == "accratio") {
            total += "| " + cell_padding(std::to_string(totals[i]/nmove), "r", max_size_col[i]) + " ";
        }
        else {
            total += "| " + cell_padding(std::to_string(static_cast<int>(totals[i])), "r", max_size_col[i]) + " ";
        }
    }
    total += "|\n";

    if (style == "basic") {
        table += std::string(head_row.size(), '-') + "\n"; 
        table += head_row + "\n";
        table += std::string(head_row.size(), '-') + "\n";
        table += rows;
        table += std::string(head_row.size(), '-') + "\n"; 
        table += total;
        table += std::string(head_row.size(), '-') + "\n"; 
    }
    else if (style == "latex") {
        table += "\begin{center}\n\begin{tabular}{";
    }
    else {
        exit(0);
    }

    if (print) {
        std::cout << std::endl;
        std::cout << table << std::endl;
        std::cout << std::endl;
    }

    return table;
}


/* ----------------------------------------------------------------------------
   Print constraint statistics
---------------------------------------------------------------------------- */

std::string System::print_constraint_statistics(const std::string &style, bool print)
{
    std::vector<std::string> keywords = {"Constraint", "#reject", "CPU-time (s)", "box-ID"};
    std::vector<std::vector<std::string> > values;
    for (Box *box : boxes) {
        for (Constraint *constraint : box->constraints) {
            values.push_back({});
            std::size_t vs = values.size();
            std::cout << constraint->label << std::endl;
            values[vs].push_back(constraint->label);
            values[vs].push_back(std::to_string(constraint->nreject));
            values[vs].push_back(std::to_string(constraint->cum_time));
            values[vs].push_back(std::to_string(box->box_id));
            std::cout << values[vs].size() << std::endl;
        }
    }
    return create_table(keywords, values, style, print);
}


/* ----------------------------------------------------------------------------
   System destructor, free memory allocated within this class
---------------------------------------------------------------------------- */

System::~System()
{
    if (!rng_allocated_externally) {
        delete rng;
    }
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    for (unsigned int i=nbox; i--;) {
        rm_box(i);
    }
    for (unsigned int i=nmove; i--;) {
        rm_move(i);
    }
    //chdir(original_dir);
}
