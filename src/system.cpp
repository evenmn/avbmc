#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

#include <mpi.h>

#include "system.h"
#include "box.h"
#include "tqdm.h"
#include "dump.h"
#include "thermo.h"
#include "rng/rng.h"
#include "boundary/boundary.h"
#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
//#include "integrator/velocityverlet.h"
#include "sampler/metropolis.h"
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

System::System(std::string working_dir_in)
{
    working_dir = working_dir_in;
    time = temp = chempot = 0.;

    logo_printed = false;
    nbox = nmove = nprocess = step = 0;
    ndim = 3;

    // initialize MPI
    int initialized_mpi;
    MPI_Initialized(&initialized_mpi);
    if (!initialized_mpi) {
        MPI_Init(nullptr, nullptr);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

void System::set_mass(const std::string label, const double mass)
{
    mass_labels.push_back(label);
    masses.push_back(mass);
}


/* ----------------------------------------------------------------------------
   Overwrite default integrator object, velocity Verlet
------------------------------------------------------------------------------- */
/*
void System::set_integrator(class Integrator* integrator_in)
{
    integrator = integrator_in;
}
*/

/* ----------------------------------------------------------------------------
   Overwrite default sampler, Metropolis
------------------------------------------------------------------------------- */

void System::set_sampler(class Sampler* sampler_in)
{
    sampler = sampler_in;
}


/* ----------------------------------------------------------------------------
   Overwrite default random number generator, Mersenne Twister
------------------------------------------------------------------------------- */

void System::set_rng(class RandomNumberGenerator* rng_in)
{
    rng = rng_in;
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

void System::add_box(Box* box_in)
{
    box_in->box_id = nbox;
    boxes.push_back(box_in);
    nbox ++;
}


/* ----------------------------------------------------------------------------
   The particles in the system have to be a subset of the particles that are
   given mass and type. That has to be checked after parsing, but before 
   simulation is started.
------------------------------------------------------------------------------- */

void System::check_masses()
{
    //for (std::string mass_label : mass_labels) {
    //    forcefield->label2type.at(mass_label);
    //}
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
    maxiter = step + nsteps/nprocess + nsteps % nprocess;
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

    if (rank == 0) {
        if (!logo_printed) { 
            print_logo();
            print_info();
        }
        print_mc_info();
        std::cout << std::endl;
        std::cout << "            Running Monte Carlo Simulation" << std::endl;
        std::cout << "=======================================================" << std::endl;
    }

    int maxiter = get_maxiter(nsteps);

    // run Monte Carlo simulation
    double start = MPI_Wtime();
    tqdm bar;
    while (step < maxiter) {
        if (rank == 0){
            bar.progress(step * nprocess, maxiter * nprocess);
        }
        for (Box* box : boxes) {
            box->dump->print_frame(step);
            box->thermo->print_line(step);
        }
        sampler->sample(nmoves);
        step ++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

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
}


/* ----------------------------------------------------------------------------
   System destructor, finalizing MPI
------------------------------------------------------------------------------- */

System::~System()
{
    MPI_Finalize();
}
