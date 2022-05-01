#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>
#include <functional>

//#include <mpi.h>

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
#include "sampler/metropolis.h"
#include "sampler/umbrella.h"
#include "moves/moves.h"
#include "particle.h"
#include "distance_manager.h"


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
    nbox = nmove = step = rank = 0;
    time = temp = chempot = poteng = 0.;
    nprocess = 1;
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

    /*
    // initialize MPI
    int initialized_mpi;
    MPI_Initialized(&initialized_mpi);
    if (!initialized_mpi) {
        MPI_Init(nullptr, nullptr);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    */
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
    rank = other.rank;
    //nprocess = other.nprocess;
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
                }
            }
            else {
                if (boxes[box_id]->forcefield_allocated_in_system) {
                    delete boxes[box_id]->forcefield;
                }
                boxes[box_id]->forcefield = new IdealGas(boxes[box_id], labels);
                boxes[box_id]->forcefield_allocated_in_system = true;
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
                }
            }
            else {
                if (boxes[box_id]->forcefield_allocated_in_system) {
                    delete boxes[box_id]->forcefield;
                }
                boxes[box_id]->forcefield = new LennardJones(boxes[box_id], paramfile);
                boxes[box_id]->forcefield_allocated_in_system = true;
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
                }
            }
            else {
                if (boxes[box_id]->forcefield_allocated_in_system) {
                    delete boxes[box_id]->forcefield;
                }
                boxes[box_id]->forcefield = new Vashishta(boxes[box_id], paramfile);
                boxes[box_id]->forcefield_allocated_in_system = true;
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

void System::add_constraint(class Constraint* constraint_in, int box_id)
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
    //    MPI_Abort(MPI_COMM_WORLD, 143);
    //}
    nmove ++;
    moves.push_back(move);
    moves_prob.push_back(prob);
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
   Returns the last iteration
------------------------------------------------------------------------------- */

int System::get_maxiter(const int nsteps)
{
    int maxiter = 0;
    if (step == 0) {
        maxiter += 1;
    }
    maxiter += step + nsteps/nprocess + nsteps % nprocess;
    return maxiter;
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
    std::cout << "Started computation on " << std::ctime(&start_time)
              << "Running on " << nprocess << " CPU threads using MPI" << std::endl;

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
   Run Monte Carlo simulation
------------------------------------------------------------------------------- */

void System::run_mc(const int nsteps, const int nmoves)
{
    for (Box* box : boxes) {
        box->nsystemsize.resize(box->npar + 1);
        box->nsystemsize[box->npar] ++;
        std::cout << box->store_distance << std::endl;
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

    /*
    if (rank == 0) {
        if (!logo_printed) { 
            print_info();
        }
        print_mc_info();
        std::cout << std::endl;
        std::cout << "            Running Monte Carlo Simulation" << std::endl;
        std::cout << "=======================================================" << std::endl;
    }
    */

    int maxiter = get_maxiter(nsteps);

    // run Monte Carlo simulation
    //double start = MPI_Wtime();
    tqdm bar;
    while (step < maxiter) {
        //if (rank == 0){
        bar.progress(step * nprocess, maxiter * nprocess);
        //}
        for (Box* box : boxes) {
            box->dump->print_frame(step);
            box->thermo->print_line(step);
        }
        sampler->sample(nmoves);
        step ++;
    }
    std::cout << std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //double end = MPI_Wtime();

    // write acceptance information to terminal
    /*
    if (rank == 0) {
        std::cout << "\n" << std::endl;
        std::cout << "Time elapsed during the job: " << end - start << "s" << std::endl;
        std::cout << std::endl;
        std::cout << "                      Acceptance Ratio" << std::endl;
        std::cout << "===================================================================================" << std::endl;
        std::cout << "Move type\t #attempts\t #accepted\t acceptance ratio\t time (s)" << std::endl;
    }
    for(Moves* move : moves){
        int ndrawntot, naccepttot;
        MPI_Reduce(&move->ndrawn, &ndrawntot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&move->naccept, &naccepttot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << move->label << "\t "
                      << ndrawntot << "\t\t "
                      << naccepttot << "\t\t "
                      << (double) naccepttot / ndrawntot << "\t\t "
                      << move->cum_time
                      << std::endl;
        }
    }
    if (rank == 0) {
        std::cout << "===================================================================================" << std::endl;
    }
    */
}


/* ----------------------------------------------------------------------------
   System destructor, finalizing MPI
------------------------------------------------------------------------------- */

System::~System()
{
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
    }
    //MPI_Finalize();
}
