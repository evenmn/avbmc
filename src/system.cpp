#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>
#include <memory>

//#include <mpi.h>

#include "system.h"
#include "box.h"
#include "tqdm.h"
#include "boundary/boundary.h"
#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
#include "rng/mersennetwister.h"
//#include "integrator/velocityverlet.h"
#include "sampler/metropolis.h"
#include "moves/moves.h"
#include "particle.h"
#include "molecule.h"


/* -------------------------------------------------------
   System constructor, taking the working directory
   'working_dir_in' as argument
---------------------------------------------------------- */

System::System(std::string working_dir_in)
{
    working_dir = working_dir_in;
    time = temp = chempot = 0.;

    initialized = false;
    nbox = ntype = nmove = nprocess = step = 0;
    ndim = 3;

    // set default objects
    //MersenneTwister rngrng();
    //rng = &rngrng;
    //integrator = new VelocityVerlet(this);
    //forcefield = new LennardJones(this);
    //sampler = new Metropolis(this);
    MoleculeTypes moleculetypes(this);
    molecule_types = &moleculetypes;
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


/* --------------------------------------------------------
   Set box temperature
----------------------------------------------------------- */

void System::set_temp(const double temp_in)
{
    temp = temp_in;
}


/* --------------------------------------------------------
   Set chemical potential of system
----------------------------------------------------------- */

void System::set_chempot(const double chempot_in)
{
    chempot = chempot_in;
}


/* --------------------------------------------------------
   Set mass of chemical symbol. Masses of all chemical symbols
   have to be given if running molecular dynamics simulations,
   as the software does not look up the masses in a table.
----------------------------------------------------------- */

void System::set_mass(const std::string label, const double mass)
{
    mass_labels.push_back(label);
    masses.push_back(mass);
}


/* --------------------------------------------------------
   Overwrite default forcefield object
----------------------------------------------------------- */

void System::set_forcefield(class ForceField* forcefield_in)
{
    forcefield = forcefield_in;
}


/* --------------------------------------------------------
   Overwrite default forcefield object, velocity Verlet
----------------------------------------------------------- */
/*
void System::set_integrator(class Integrator* integrator_in)
{
    integrator = integrator_in;
}
*/

/* --------------------------------------------------------
   Overwrite default sampler, Metropolis
----------------------------------------------------------- */

void System::set_sampler(class Sampler* sampler_in)
{
    sampler = sampler_in;
}


/* --------------------------------------------------------
   Overwrite default random number generator, Mersenne
   Twister
----------------------------------------------------------- */

void System::set_rng(class RandomNumberGenerator* rng_in)
{
    rng = rng_in;
}


/* --------------------------------------------------
   Add move type and the corresponding probability.
   The probabilities have to add up to 1.
----------------------------------------------------- */
/*
void System::add_move(Moves move, double prob)
{
    nmove ++;
    moves.emplace_back(move);
    moves_prob.push_back(prob);
}
*/
/*
void System::add_move(Moves& move, double prob)
{
    nmove ++;
    moves.emplace_back((*move));
    moves_prob.push_back(prob);
}
*/

void System::add_move(Moves* move, double prob)
{
    nmove ++;
    moves.emplace_back(move);
    moves_prob.push_back(prob);
}


/* ------------------------------------------------------
   Add new molecule type when molecule is defined by
   one atom.
--------------------------------------------------------- */

void System::add_molecule_type(std::string element, const double molecule_prob)
{
    std::vector<std::string> elements = {element};
    std::valarray<double> default_atom(0.0, ndim);
    std::vector<std::valarray<double> > default_mol = {default_atom};
    molecule_types->add_molecule_type(elements, 0.0, molecule_prob, default_mol);
}


/* ------------------------------------------------------
   Add new molecule type consisting of several atoms
--------------------------------------------------------- */

void System::add_molecule_type(std::vector<std::string> elements, const double rc, const double molecule_prob,
                            std::vector<std::valarray<double> > default_mol)
{
    molecule_types->add_molecule_type(elements, rc, molecule_prob, default_mol);
}


/* -----------------------------------------------------
   Add box 'box_in' to system
-------------------------------------------------------- */

void System::add_box(Box* box_in)
{
    nbox ++;
    boxes.push_back(box_in);
}
/*
void System::add_box(std::shared_ptr<Box> box_in)
{
    boxes.emplace_back(box_in);
    nbox ++;
}

void System::add_box(Box* box_in)
{
    boxes.emplace_back(box_in);
    nbox ++;
}
*/

/* -----------------------------------------------------
   The particles in the system have to be a subset of 
   the particles that are given mass and  type. That 
   has to be checked after parsing, but before 
   simulation is started.
-------------------------------------------------------- */

void System::check_masses()
{
    // Check that all particles are assigned a mass
    for (std::string mass_label : mass_labels){
        label2type.at(mass_label);
    }
    /*
    for (std::string unique_label : unique_labels){
        bool label_covered = false;
        for (std::string mass_label : mass_labels){
            if (unique_label == mass_label){
                label_covered = true;
            }
        }
        //std::string msg = "Particle type " + unique_label + " is not assigned a mass! Aborting.";
        assert (label_covered);
    }
    */
}


/* ----------------------------------------------------
   Initialize molecules used by AVBMCMol types of moves
------------------------------------------------------- */

void System::init_molecules()
{
    // If molecule configuration is not set, let all single particles
    // be considered as molecules
    if (!molecule_types->configured)
    {
        for(std::string label : unique_labels){
            add_molecule_type(label, 1./ntype);
        }
    }

    // Since particles are associated with types rather than labels,
    // also molecule configuration labels have to be converted to types
    for(std::vector<std::string> elements : molecule_types->molecule_elements){
        std::vector<int> types;
        for (int i=0; i < ntype; i++){
            for(std::string element : elements){
                if(element == unique_labels[i]){
                    types.push_back(i);
                }
            }
        }
        molecule_types->molecule_types.push_back(types); 
    }
}


/* -------------------------------------------------------
   Initialize variables needed before simulation.
   Ensure that all particles are covered by parameter file.
   All forcefields should have a label1_vec, which covers
   all element labels
---------------------------------------------------------- */

void System::init_simulation()
{
    // find unique labels
    unique_labels = forcefield->label1_vec;
    std::sort( unique_labels.begin(), unique_labels.end() );
    unique_labels.erase( std::unique( unique_labels.begin(),
                         unique_labels.end() ), unique_labels.end() );
    ntype = unique_labels.size();

    // map labels to types
    for (int i=0; i < ntype; i++){
        label2type[unique_labels[i]] = i;
    }

    // go through all particles in all boxes and assert that    
    // all element labels are covered by parameter file
    for (Box* box : boxes){
        for (Particle& particle : box->particles){
            try {
                // will throw exception if label is not known
                particle.type = label2type.at(particle.label);
                std::cout << "PARTICLE TYPE: " << particle.type << std::endl;
            }
            catch (std::out_of_range) {
                std::cout << "Particle label '" + particle.label + "' is not covered by parameter file" << std::endl;
                //MPI_Abort(MPI_COMM_WORLD, 143);
            }
            /*
            bool particle_covered = false;
            for (int j=0; j < ntype; j++) {
                if (particle->label == unique_labels[j]) {
                    particle->type = j;
                    particle_covered = true;
                }
            }
            //std::string msg = "Particle type " + particle->label + " is not covered by parameter file! Aborting.";
            assert (particle_covered);
            */
        }
    }
    for (Box* box : boxes) {
        for (Particle particle : box->particles) {
            std::cout << "PARTICLE TYPES 2: " << particle.type << std::endl;
        }
    }
    // Sort forcefield parameters according to particle types
    forcefield->sort_params();
}


/* -------------------------------------------------------
   Returns the last iteration
---------------------------------------------------------- */

int System::get_maxiter(const int nsteps)
{
    //int maxiter;
    int maxiter = step + nsteps/nprocess + nsteps % nprocess;
    //if(step == 0){
    //    maxiter = nsteps + 1;
    //}
    //else{
    //    maxiter = nsteps + step;
    //}
    return maxiter;
}


/* -------------------------------------------------------
   Print logo header
---------------------------------------------------------- */

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
}


/* -------------------------------------------------------
   Print information about simulation
---------------------------------------------------------- */

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
    std::cout << "Forcefield:               " << forcefield->label 
              << " " << forcefield->paramfile << std::endl;
    std::cout << "Random number generator:  " << rng->label << std::endl;
    std::cout << std::endl;
    std::cout << "Number of particle types: " << ntype << std::endl;
    std::cout << "Unique particle types:";
    for(int i=0; i < ntype; i++){
        std::cout << " " << unique_labels[i];
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Number of boxes:          " << nbox << std::endl;
    for(int i=0; i < nbox; i++){
        std::cout << "  Box " << i+1 << ":" << std::endl;
        std::cout << "    Number of atoms:   " << boxes[i]->npar << std::endl;
        std::cout << "    Boundary:          " << boxes[i]->boundary->label << std::endl;
    }
    std::cout << std::endl;
    /*
    std::cout << "Number of molecule types: " << molecule_types->ntype << std::endl;
    for(int i=0; i < molecule_types->ntype; i++){
        std::cout << "  Molecule type " << i+1 << ":" << std::endl;
        std::cout << "    Atoms:            ";
        for(std::string element : molecule_types->molecule_elements[i]){
            std::cout << " " << element;
        }
        std::cout << std::endl;
        std::cout << "    Critical distance: " << molecule_types->rcs[i] << std::endl;
        std::cout << "    Probability:       " << molecule_types->molecule_probs[i] << std::endl;
    }
    */
    std::cout << std::endl;
}


/* --------------------------------------------------------
   Print Monte Carlo information
----------------------------------------------------------- */

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


/*
void Box::run_md(const int nsteps)
{
    // Run molecular dynamics simulation
    //
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
        integrator->next_step();
        step ++;
    }
    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = chrono::duration<double>(t2 - t1).count();
    cout << "Elapsed time: " << duration_seconds << "s" << endl;
}
*/

/* -------------------------------------------------------
   Run Monte Carlo simulation
---------------------------------------------------------- */

void System::run_mc(const int nsteps, const int nmoves)
{
    init_simulation();
    //init_molecules();
    for(Box* box : boxes){
        box->nsystemsize.resize(box->npar + 1);
        box->nsystemsize[box->npar] ++;
    }

    for (Moves* move : moves){
        move->ndrawn = 0;
        move->naccept = 0;
    }

    if (rank == 0) {
        print_logo();
        print_info();
        print_mc_info();
        std::cout << std::endl;
        std::cout << "            Running Monte Carlo Simulation" << std::endl;
        std::cout << "=======================================================" << std::endl;
    }
    //thermo->print_header();

    int maxiter = nsteps; //get_maxiter(nsteps);
    std::cout << "After get_maxiter" << std::endl;

    // run Monte Carlo simulation
    //double start = MPI_Wtime();
    tqdm bar;
    //bar.set_theme_circle();
    while(step < maxiter){
        std::cout << step << std::endl;
        //if (rank == 0){
        //    bar.progress(step * nprocess, maxiter * nprocess);
        //}
        //box->dump->print_frame(step);
        //box->thermo->print_line(step);
        sampler->sample(nmoves);
        step ++;
    }
    /*
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    if (rank == 0) {
        std::cout << "\n" << std::endl;
        std::cout << "Time elapsed during the job: " << end - start << "s" << std::endl;
        std::cout << std::endl;
        std::cout << "                      Acceptance Ratio" << std::endl;
        std::cout << "=================================================================" << std::endl;
        std::cout << "Move type\t #drawn\t\t #accepted\t acceptance ratio" << std::endl;
    }
    for(Moves move : moves){
        int ndrawntot, naccepttot;
        MPI_Reduce(&move.ndrawn, &ndrawntot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&move.naccept, &naccepttot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << move.label << "\t "
                      << ndrawntot << "\t\t "
                      << naccepttot << "\t\t "
                      << (double) naccepttot / ndrawntot
                      << std::endl;
        }
    }
    if (rank == 0) {
        std::cout << "=================================================================" << std::endl;
    }
    */
}


/* ---------------------------------------------------------
   System destructor, finalizing MPI
------------------------------------------------------------ */
System::~System()
{
    //delete rng;
    //delete forcefield;
    //delete sampler;
    //delete molecule_types;
    //rng = nullptr;
    //forcefield = nullptr;
    //sampler = nullptr;
    //molecule_types = nullptr;
    //MPI_Finalize();
}
