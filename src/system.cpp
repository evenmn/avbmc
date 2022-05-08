#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>
#include <functional>

#include "system.h"
#include "box.h"
#include "tqdm.h"
#include "dump.h"
#include "thermo.h"
#include "rng/rng.h"
#include "rng/mersennetwister.h"
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
#include "particle.h"
#include "distance_manager.h"
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
------------------------------------------------------------------------------- */

System::System(const std::string &working_dir_in, bool initialize_in) 
    : working_dir(working_dir_in)
{
    nbox = nmove = step = 0;
    time = temp = chempot = poteng = 0.;
    ndim = 3;

    //print_logo();
    logo_printed = rng_allocated_externally = sampler_allocated_externally = false;

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
------------------------------------------------------------------------------- */

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
   Set box temperature used in NVT and uPT Monte Carlo simulations
------------------------------------------------------------------------------- */

void System::set_temp(const double temp_in)
{
    temp = temp_in;
}


/* ----------------------------------------------------------------------------
   Set chemical potential of system used in grand canonical ensemble
------------------------------------------------------------------------------- */

void System::set_chempot(const double chempot_in)
{
    chempot = chempot_in;
}


/* ----------------------------------------------------------------------------
   Set seed to be used by RNG
------------------------------------------------------------------------------- */

void System::set_seed(unsigned int seed_)
{
    rng->set_seed(seed_);
}


/* ----------------------------------------------------------------------------
   Set mass of chemical symbol. Masses of all chemical symbols have to be given
   if running molecular dynamics simulations, as the software does not look up
   the masses in a table.
------------------------------------------------------------------------------- */
/*
void System::set_mass(const std::string label, const double mass)
{
    mass_labels.push_back(label);
    masses.push_back(mass);
}
*/


/* ----------------------------------------------------------------------------
   Overwrite default sampler, Metropolis
------------------------------------------------------------------------------- */

void System::set_sampler(class Sampler *sampler_in)
{
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    sampler = sampler_in;
    sampler_allocated_externally = true;
}


void System::set_sampler(const std::string &sampler_in, std::function<double(int)> f)
{
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    if (sampler_in == "metropolis") {
        sampler = new Metropolis(this);
        sampler_allocated_externally = false;
    }
    else if (sampler_in == "umbrella") {
        sampler = new Umbrella(this, f);
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
------------------------------------------------------------------------------- */

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
    if (rng_in == "mersennetwister") {
        rng = new MersenneTwister;
        rng_allocated_externally = false;
    }
    else if (rng_in == "mt19937") {
        rng = new MersenneTwister;
        rng_allocated_externally = false;
    }
    else {
        std::cout << "Random number generator '" << rng_in << "' is not implemented!"
                  << "Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Set forcefield. This method is forward to the box and can only be done if a
   box was detected.
------------------------------------------------------------------------------- */

void System::set_forcefield(class ForceField* forcefield_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! Aborting." << std::endl;
        exit(0);
    }
    else {
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
}


void System::set_forcefield(const std::string &forcefield_in,
    const std::vector<std::string> &labels, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! Aborting." << std::endl;
        exit(0);
    }
    else {
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
            std::cout << "Force-field '" << forcefield_in << "' is not implemented or"
                      << "does not have the given signature! Aborting." << std::endl;
            exit(0);
        }
    }
}


void System::set_forcefield(const std::string &forcefield_in,
    const std::string &paramfile, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set forcefield! Aborting." << std::endl;
        exit(0);
    }
    else {
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
}
   

/* ----------------------------------------------------------------------------
   Set boundary. This method is forward to the box and can only be done if a
   box was detected.
------------------------------------------------------------------------------- */

void System::set_boundary(class Boundary *boundary_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set boundary" << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->set_boundary(boundary_in);
                box->boundary_allocated_in_system = false;
            }
        }
        else {
            boxes[box_id]->set_boundary(boundary_in);
            boxes[box_id]->boundary_allocated_in_system = false;
        }
    }
}


void System::set_boundary(const std::string &boundary_in,
    std::valarray<double> length, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set boundary" << std::endl;
        exit(0);
    }
    else {
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
}


/* ----------------------------------------------------------------------------
   Add constraint. This method is forward to the box and can only be done if a
   box was detected.
------------------------------------------------------------------------------- */

void System::add_constraint(Constraint* constraint_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add constraint" << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            if (nbox > 1) {
                std::cout << "Warning: More than one box was detected. Setting"
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
}


void System::add_constraint(const std::string &constraint_in,
    const std::string &element1, const std::string &element2, double distance,
    int nneigh, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add constraint" << std::endl;
        exit(0);
    }
    else if (box_id < 0 && nbox > 1) {
        std::cout << "Warning: More than one box was detected. Setting"
                  << "constraint for all boxes" << std::endl;
    }
    else {
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
}


/* ----------------------------------------------------------------------------
   Take a snapshot of a given box.
------------------------------------------------------------------------------- */

void System::snapshot(const std::string &filename, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot take snapshot! Aborting." << std::endl;
        exit(0);
    }
    else {
        boxes[box_id]->snapshot(filename);
    }
}


/* ----------------------------------------------------------------------------
   Set dump. This method is forward to the box and can only be done if
   a box was detected.
------------------------------------------------------------------------------- */

void System::set_dump(int freq, const std::string &filename,
    const std::vector<std::string> &outputs, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set dump! Aborting." << std::endl;
        exit(0);
    }
    else {
        boxes[box_id]->set_dump(freq, filename, outputs);
    }
}


/* ----------------------------------------------------------------------------
   Set thermo. This method is forward to the box and can only be done if a box
   was detected.
------------------------------------------------------------------------------- */

void System::set_thermo(int freq, const std::string &filename,
    const std::vector<std::string> &outputs, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot set thermo! Aborting." << std::endl;
        exit(0);
    }
    else {
        boxes[box_id]->set_thermo(freq, filename, outputs);
    }
}


/* ----------------------------------------------------------------------------
   Add move type and the corresponding probability. The probabilities have to
   add up to 1. Forcefield has to initialized first, to link AVBMC atoms to
   types.
------------------------------------------------------------------------------- */

void System::add_move(Moves* move, double prob)
{
    //if (!initialized) {
    //    std::cout << "Forcefield needs to be initialized before adding moves!" << std::endl;
    //    exit(0);
    //}
    nmove ++;
    moves.push_back(move);
    moves_prob.push_back(prob);
    moves_allocated_in_system.push_back(false);
}


void System::add_move(const std::string &move_in, double prob, double dx, double Ddt, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (move_in == "trans") {
            if (box_id < 0) {
                for (Box *box : boxes) {
                    Moves *move = new Trans(this, box, dx);
                    add_move(move, prob);
                    moves_allocated_in_system[nmove-1] = true;
                }
            }
            else {
                Moves *move = new Trans(this, boxes[box_id], dx);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
        }
        else if (move_in == "transmh") {
            if (box_id < 0) {
                for (Box *box : boxes) {
                    Moves *move = new TransMH(this, box, dx, Ddt);
                    add_move(move, prob);
                    moves_allocated_in_system[nmove-1] = true;
                }
            }
            else {
                Moves *move = new TransMH(this, boxes[box_id], dx, Ddt);
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
}


void System::add_move(const std::string &move_in, double prob, 
    const std::string &particle_in, double r_below, double r_above, bool energy_bias, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (move_in == "avbmc") {
            if (box_id < 0) {
                for (Box *box : boxes) {
                    Moves *move = new AVBMC(this, box, particle_in, r_below, r_above, energy_bias);
                    add_move(move, prob);
                    moves_allocated_in_system[nmove-1] = true;
                }
            }
            else {
                Moves *move = new AVBMC(this, boxes[box_id], particle_in, r_below, r_above, energy_bias);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
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
        else {
            std::cout << "Move '" << move_in << "' is not implemented!"
                      << "Aborting." << std::endl;
            exit(0);
        }
    }
}


void System::add_move(const std::string &move_in, double prob,
    std::vector<Particle> molecule_in, double r_below,
    double r_above, double r_inner, bool energy_bias, bool target_mol, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add move! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (move_in == "avbmcmol") {
            if (box_id < 0) {
                for (Box *box : boxes) {
                    Moves *move = new AVBMCMol(this, box, molecule_in, r_below, r_above, r_inner, energy_bias, target_mol);
                    add_move(move, prob);
                    moves_allocated_in_system[nmove-1] = true;
                }
            }
            else {
                Moves *move = new AVBMCMol(this, boxes[box_id], molecule_in, r_below, r_above, r_inner, energy_bias, target_mol);
                add_move(move, prob);
                moves_allocated_in_system[nmove-1] = true;
            }
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
        else {
            std::cout << "Move '" << move_in << "' is not implemented!"
                      << "Aborting." << std::endl;
            exit(0);
        }
    }
}


/* ----------------------------------------------------------------------------
   Add box 'box_in' to system
------------------------------------------------------------------------------- */

void System::add_box()
{
    Box *box = new Box(this);
    add_box(box);
    box->box_allocated_in_system = true;
}


void System::add_box(Box* box_in)
{
    box_in->box_id = nbox;
    boxes.push_back(box_in);
    nbox ++;
}


/* ----------------------------------------------------------------------------
   Add particles to the system. 
------------------------------------------------------------------------------- */

void System::add_particle(Particle particle_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->add_particle(particle_in);
            }
        }
        else {
            boxes[box_id]->add_particle(particle_in);
        }
    }
}


void System::add_particle(const std::string &label_in, std::valarray<double> pos, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->add_particle(label_in, pos);
            }
        }
        else {
            boxes[box_id]->add_particle(label_in, pos);
        }
    }
}


void System::add_particles(std::vector<Particle> particles_in, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->add_particles(particles_in);
            }
        }
        else {
            boxes[box_id]->add_particles(particles_in);
        }
    }
}


void System::add_particles(const std::string &label_in, std::vector<std::valarray<double> > pos, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->add_particles(label_in, pos);
            }
        }
        else {
            boxes[box_id]->add_particles(label_in, pos);
        }
    }
}


void System::read_particles(const std::string &filename, int box_id)
{
    if (nbox < 1) {
        std::cout << "No box found, cannot add particles! Aborting." << std::endl;
        exit(0);
    }
    else {
        if (box_id < 0) {
            for (Box *box : boxes) {
                box->read_particles(filename);
            }
        }
        else {
            boxes[box_id]->read_particles(filename);
        }
    }
}


/* ----------------------------------------------------------------------------
   Write histogram of system sizes out to file
------------------------------------------------------------------------------- */

void System::write_size_histogram(const std::string &filename, int box_id)
{
    boxes[box_id]->write_size_histogram(filename);
}


std::vector<unsigned int> System::get_size_histogram(int box_id)
{
    return boxes[box_id]->size_histogram;
}

/* ----------------------------------------------------------------------------
   Returns the last iteration
------------------------------------------------------------------------------- */

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
------------------------------------------------------------------------------- */

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
------------------------------------------------------------------------------- */

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
    for(int i=0; i < nbox; i++){
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
    for(int i=0; i < nmove; i++){
        std::cout << "  Move type " << i+1 << ": ";
        std::cout << moves[i]->repr();
        std::cout << "    Probability: " << moves_prob[i] << std::endl;
    }
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
   Run molecular dynamics simulation
------------------------------------------------------------------------------- */
/*
void Box::run_md(const int nsteps)
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
------------------------------------------------------------------------------- */

void System::run_mc(const unsigned int nsteps, const unsigned int nmoves)
{
    for (Box* box : boxes) {
        box->size_histogram.resize(box->npar + 1);
        box->size_histogram[box->npar] ++;
        if (box->store_distance) {
            box->distance_manager->initialize();
        }
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


    unsigned int init_step = step;
    tqdm bar;
    unsigned int maxiter = get_maxiter(nsteps);
    while (step < maxiter) {
        bar.progress(step-init_step, nsteps);
        for (Box* box : boxes) {
            box->dump->print_frame(step);
            box->thermo->print_line(step);
        }
        sampler->sample(nmoves);
        step ++;
    }
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
   System destructor, free memory allocated within this class
------------------------------------------------------------------------------- */

System::~System()
{
    unsigned int i;

    if (!rng_allocated_externally) {
        delete rng;
    }
    if (!sampler_allocated_externally) {
        delete sampler;
    }
    for (Box *box : boxes) {
        if (box->forcefield_allocated_in_system) {
            delete box->forcefield;
        }
        if (box->boundary_allocated_in_system) {
            delete box->boundary;
        }
        if (box->box_allocated_in_system) {
            delete box;
        }
        for (i=0; i<box->nconstraint; i++) {
            if (box->constraint_allocated_in_system[i]) {
                delete box->constraints[i];
            }
        }
    }
    for (i=0; i<nmove; i++) {
        if (moves_allocated_in_system[i]) {
            delete moves[i];
        }
    }
}
