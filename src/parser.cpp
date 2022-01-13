#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <functional> 
#include <cassert>
#include <iterator>
#include <vector>

#include "parser.h"
#include "io.h"
#include "box.h"
#include "system.h"
#include "init_position.h"
#include "particle.h"
#include "rng/mersennetwister.h"
//#include "init_velocity.h"

//#include "integrator/euler.h"
//#include "integrator/eulercromer.h"
//#include "integrator/velocityverlet.h"
//#include "integrator/rungekutta4.h"

#include "forcefield/lennardjones.h"
#include "forcefield/vashishta.h"

#include "sampler/metropolis.h"
#include "sampler/umbrella.h"

#include "moves/trans.h"
#include "moves/transmh.h"
#include "moves/avbmc.h"
#include "moves/avbmcin.h"
#include "moves/avbmcout.h"
#include "moves/avbmcmol.h"
#include "moves/avbmcmolin.h"
#include "moves/avbmcmolout.h"

#include "boundary/fixed.h"
#include "boundary/stillinger.h"



/* -----------------------------------------------
   Primary parser. Checks what kind of command
   that is given (set, add, take, run, thermo
   or dump) and calls the correct function
-------------------------------------------------- */

void parser(int argc, char** argv)
{
    if (argc >= 2) {  // check if input script is given
        std::string filename = argv[1];
        std::ifstream f(filename);

        if (f.is_open()) { 
            System system("");
            Box box(&system);
            Metropolis sampler(&system);
            MersenneTwister rng(void);
            system.add_box(&box);
            system.set_sampler(&sampler);
            system.set_rng(&rng);
            std::string line;
            while (std::getline(f, line))
            {
                if (line.empty()) {
                    // Continue if line is blank
                    continue;
                }
                else if (line.rfind("#", 0) == 0) {
                    // Continue if line is commented out
                    continue;
                }
                else {
                    std::vector<std::string> splitted = split(line);
                    int argc = splitted.size();
                    std::string cmd_cat = splitted[0];
                    if (cmd_cat == "set") {
                        set(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "add") {
                        add(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "take") {
                        take(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "run") {
                        run(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "thermo") {
                        thermo(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "dump") {
                        dump(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "write") {
                        write(system, box, splitted, argc);
                    }
                    else if (cmd_cat == "rm") {
                        rm(system, box, splitted, argc);
                    }
                    else {
                        std::cout << "There is no category '" + cmd_cat + "'! Aborting." << std::endl;
                        exit(0);
                    }
                }
            }
        }
        else {
            std::cout << "Could not open the file '" + filename + "'!" << std::endl;
            exit(0);
        }
    }
    else {  // printing some information if no input script is given
        std::cout << "avbmc version 0.0.1" << std::endl;
        std::cout << "Author: Even M. Nordhagen" << std::endl;
        std::cout << "url: github.com/evenmn/avbmc" << std::endl;
    }
}


/* ----------------------------------------------------
   This function takes care of the set-commands. The 
   set commands are used to set box properties, and
   overwrite the previous value
------------------------------------------------------- */

void set(System & system, Box& box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 2);
    std::string keyword = splitted[1];
    if(keyword == "temp"){
        system.set_temp(std::stod(splitted[2]));
    }
    else if(keyword == "chempot"){
        system.set_chempot(std::stod(splitted[2]));
    }
    else if(keyword == "mass"){
        assert(argc > 3);
        std::string label = splitted[2];
        double mass = std::stod(splitted[3]);
        system.set_mass(label, mass);
    }
    else if(keyword == "forcefield"){
        assert(argc > 3);
        std::string forcefield = splitted[2];
        std::string paramfile = splitted[3];
        if(forcefield == "lennardjones"){
            LennardJones forcefield(&system, paramfile);
            system.set_forcefield(&forcefield);
        }
        else if(forcefield == "vashishta"){
            Vashishta forcefield(&system, paramfile);
            system.set_forcefield(&forcefield);
        }
        else{
            std::cout << "Forcefield '" + forcefield + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    /*
    else if(keyword == "integrator"){
        std::string integrator = splitted[2];
        if(argc == 3){
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box));
            }
            else{
                std::cout << "Integrator '" + integrator + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else if(argc >= 4){
            double dt = std::stod(splitted[3]);
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box, dt));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box, dt));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box, dt));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box, dt));
            }
            else{
                std::cout << "Integrator '" + integrator + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else{
            std::cout << "Integrator does not have a sufficient number of arguments! Aborting." << std::endl;
            exit(0);
        }
    }
    */
    else if(keyword == "sampler"){
        std::string sampler = splitted[2];
        if(sampler == "metropolis"){
            Metropolis sampler(&system);
            system.set_sampler(&sampler);
        }
        else if(sampler == "umbrella"){
            assert(argc > 6);
            std::string umbrella_type = splitted[3];
            if(umbrella_type == "poly"){
                double a = std::stod(splitted[4]);
                double b = std::stod(splitted[5]);
                double c = std::stod(splitted[6]);
                auto f = [a, b, c] (const int x) { return (a * x * x + b * x + c); };
                Umbrella sampler(&system, f);
                system.set_sampler(&sampler);
            }
            else if(umbrella_type == "square"){
                double nc = std::stod(splitted[4]);
                double k = std::stod(splitted[5]);
                double nv = std::stod(splitted[6]);
                auto f = [nc, k, nv] (const int n) { return (k * (n - nc) * (n - nc)) + nv; };
                Umbrella sampler(&system, f);
                system.set_sampler(&sampler);
            }
            else{
                std::cout << "Umbrella type '" + umbrella_type + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else {
            std::cout << "Sampler '" + sampler + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "rng"){
        std::string rng = splitted[2];
        if(rng == "mersennetwister"){
            MersenneTwister rng(void);
            system.set_rng(&rng);
        }
        else{
            std::cout << "Random number generator '" + rng + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "boundary"){
        assert(argc > 3);
        std::string boundary = splitted[2];
        if(boundary == "stillinger"){
            double rc = std::stod(splitted[3]);
            Stillinger boundary(&box, rc);
            box.set_boundary(&boundary);
        }
        else if(boundary == "fixed"){
            if(argc == 4){
                // one box length, indicating one dimension
                double lx = stod(splitted[3]);
                Fixed boundary(&box, {lx});
                box.set_boundary(&boundary);
            }
            else if(argc == 5){
                // two box lengths, indicating two dimensions
                double lx = stod(splitted[3]);
                double ly = stod(splitted[4]);
                Fixed boundary(&box, {lx, ly});
                box.set_boundary(&boundary);
            }
            else if(argc == 6){
                double lx = stod(splitted[3]);
                double ly = stod(splitted[4]);
                double lz = stod(splitted[5]);
                Fixed boundary(&box, {lx, ly, lz});
                box.set_boundary(&boundary);
            }
        }
        else{
            std::cout << "Boundary '" + boundary + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    /*
    else if(keyword == "velocity"){
        assert(argc > 3);
        std::string init_style = splitted[2];
        if(init_style == "temp"){
            double temp = std::stod(splitted[3]);
            box.velocity = new Temp(box.rng, temp);
            box.velocities = box.velocity->get_velocity(box.npar, box.ndim);
        }
        else if(init_style == "gauss"){
            assert(argc > 4);
            double mean = std::stod(splitted[3]);
            double var = std::stod(splitted[4]);
            box.velocity = new Gauss(box.rng, mean, var);
            box.velocities = box.velocity->get_velocity(box.npar, box.ndim);
        }
        else{
            std::cout << "Velocity style '" + init_style + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    */
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------
   This functions takes care of the add-commands in the input
   script. The add commands are used to append a property to 
   a vector (particles, moves, molecule types etc..)
------------------------------------------------------------------- */

void add(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 1);
    std::string keyword = splitted[1];

    if(keyword == "move"){
        assert(argc > 3);
        std::string move = splitted[2];
        double prob = std::stod(splitted[3]);
        if(move == "trans"){
            if(argc == 4){
                Trans move(&system, &box);
                system.add_move(&move, prob);
            }
            else{
                double dx = std::stod(splitted[4]);
                Trans move(&system, &box, dx);
                system.add_move(&move, prob);
            }
        }
        else if(move == "transmh"){
            if(argc == 4){
                TransMH move(&system, &box);
                system.add_move(&move, prob);
            }
            else if(argc == 5){
                double dx = std::stod(splitted[4]);
                TransMH move(&system, &box, dx);
                system.add_move(&move, prob);
            }
            else{
                double dx = std::stod(splitted[4]);
                double Ddt = std::stod(splitted[5]);
                TransMH move(&system, &box, dx, Ddt);
                system.add_move(&move, prob);
            }
        }
        else if(move == "avbmc"){
            assert(argc > 4);
            std::string label = splitted[4];
            if(argc == 5){
                AVBMCIn move1(&system, &box, label); 
                AVBMCOut move2(&system, &box, label);
                system.add_move(&move1, prob / 2.);
                system.add_move(&move2, prob / 2.);
            }
            else if(argc == 6){
                double r_below = std::stod(splitted[5]);
                AVBMCIn move1(&system, &box, label, r_below);
                AVBMCOut move2(&system, &box, label);
                system.add_move(&move1, prob / 2.);
                system.add_move(&move2, prob / 2.);
            }
            else{
                double r_below = std::stod(splitted[5]);
                double r_above = std::stod(splitted[6]);
                AVBMCIn move1(&system, &box, label, r_below, r_above);
                AVBMCOut move2(&system, &box, label, r_above);
                system.add_move(&move1, prob / 2.);
                system.add_move(&move2, prob / 2.);
            }
        }
        else if(move == "avbmcin"){
            assert(argc > 4);
            std::string label = splitted[4];
            if(argc == 5){
                AVBMCIn move(&system, &box, label); 
                system.add_move(&move, prob);
            }
            else if(argc == 6){
                double r_below = std::stod(splitted[5]);
                AVBMCIn move(&system, &box, label, r_below); 
                system.add_move(&move, prob);
            }
            else{
                double r_below = std::stod(splitted[5]);
                double r_above = std::stod(splitted[6]);
                AVBMCIn move(&system, &box, label, r_below, r_above); 
                system.add_move(&move, prob);
            }
        }
        else if(move == "avbmcout"){
            assert(argc > 4);
            std::string label = splitted[4];
            if(argc == 5){
                AVBMCOut move(&system, &box, label);
                system.add_move(&move, prob);
            }
            else{
                double r_above = std::stod(splitted[5]);
                AVBMCOut move(&system, &box, label, r_above);
                system.add_move(&move, prob);
            }
        }
        /*
        else if(move == "avbmcmol"){
            if(argc == 4){
                box.add_move(new AVBMCMol(&box), prob);
            }
            else if(argc == 5){
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMCMol(&box, r_below), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                double r_above = std::stod(splitted[5]);
                box.add_move(new AVBMCMol(&box, r_below, r_above), prob);
            }
        }
        else if(move == "avbmcmolin"){
            if(argc == 4){
                box.add_move(new AVBMCInMol(&box), prob);
            }
            else if(argc == 5){
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMCInMol(&box, r_below), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                double r_above = std::stod(splitted[5]);
                box.add_move(new AVBMCInMol(&box, r_below, r_above), prob);
            }
        }
        else if(move == "avbmcmolout"){
            if(argc == 4){
                box.add_move(new AVBMCOutMol(&box), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMCOutMol(&box, r_below), prob);
            }
        }
        */
        else{
            std::cout << "Move '" + move + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "particle"){
        assert(argc > 3);
        std::string label = splitted[2];
        double x = std::stod(splitted[3]);
        if(argc == 4){
            // assuming 1d
            box.add_particle(label, {x});
        }
        if(argc == 5){
            // assuming 2d
            double y = std::stod(splitted[4]);
            box.add_particle(label, {x, y});
        }
        else{
            // assuming 3d
            double y = std::stod(splitted[4]);
            double z = std::stod(splitted[5]);
            box.add_particle(label, {x, y, z});
        }
    }
    else if(keyword == "particles"){
        assert(argc > 3);
        std::string init_type = splitted[2];
        if(init_type == "fcc"){
            assert(argc > 5);
            std::string label = splitted[3];
            int ncell = std::stoi(splitted[4]);
            double lenbulk = std::stod(splitted[5]);
            std::vector<std::valarray<double> > positions;
            if(argc == 6){
                positions = fcc(ncell, lenbulk);
            }
            else{
                int ndim = std::stoi(splitted[6]);
                positions = fcc(ncell, lenbulk, ndim);
            }
            for(std::valarray<double> position : positions){
                box.add_particle(label, position);
            }
        }
        else if(init_type == "xyz"){
            std::string filename = splitted[3];
            std::vector<Particle> particles_in = read_xyz(filename);
            box.add_particles(particles_in);
        }
        else{
            std::cout << "Particle initialization type '" + init_type + "' is not known! Aborting." << std::endl;
        }
    }
    /*
    else if (keyword == "moleculetype"){
        assert (argc > 3);
        std::vector<std::string> elements;
        std::vector<std::valarray<double> > default_mol;
        elements.push_back(splitted[2]);
        if (argc == 4){  // one atom, ex. H 1.0
            double molecule_prob = std::stod(splitted[3]);
            box.add_molecule_type(splitted[2], molecule_prob);
        }
        else if (argc == 5){ // read default molecule from xyz-file
            elements.clear();
            double rc = std::stod(splitted[2]);
            double molecule_prob = std::stod(splitted[3]);
            std::string filename = splitted[4];
            std::vector<Particle *> default_particles = read_xyz(filename);
            for(Particle* particle : default_particles){
                elements.push_back(particle->label);
                default_mol.push_back(particle->r);
            }
            box.add_molecule_type(elements, rc, molecule_prob, default_mol);
        }
        else if (argc == 8){ // two atoms, one dimension, ex. O O 1.5 1.0 0.0 1.0
            elements.push_back(splitted[3]);
            double rc = std::stod(splitted[4]);
            double molecule_prob = std::stod(splitted[5]);
            double x1 = std::stod(splitted[6]);
            double x2 = std::stod(splitted[7]);
            default_mol.push_back({x1});
            default_mol.push_back({x2});
            box.add_molecule_type(elements, rc, molecule_prob, default_mol);
        }
        else if (argc == 10){ // two atoms, two dimensions or three atoms, one dimension
            // ex. O O 1.5 1.0 0.0 0.0 1.0 0.0
            // or  O H H 1.5 1.0 -1.0 0.0 1.0
            try {
                elements.push_back(splitted[3]);
                double rc = std::stod(splitted[4]);
                double molecule_prob = std::stod(splitted[5]);
                double x1 = std::stod(splitted[6]);
                double y1 = std::stod(splitted[7]);
                double x2 = std::stod(splitted[8]);
                double y2 = std::stod(splitted[9]);
                default_mol.push_back({x1, y1});
                default_mol.push_back({x2, y2});
                box.add_molecule_type(elements, rc, molecule_prob, default_mol);
            }
            catch(...) {
                elements.push_back(splitted[3]);
                elements.push_back(splitted[4]);
                double rc = std::stod(splitted[5]);
                double molecule_prob = std::stod(splitted[6]);
                double x1 = std::stod(splitted[7]);
                double x2 = std::stod(splitted[8]);
                double x3 = std::stod(splitted[9]);
                default_mol.push_back({x1});
                default_mol.push_back({x2});
                default_mol.push_back({x3});
                box.add_molecule_type(elements, rc, molecule_prob, default_mol);
            }
        }
        else if (argc == 12){ // two atoms, three dimensions
            elements.push_back(splitted[3]);
            elements.push_back(splitted[4]);
            double rc = std::stod(splitted[4]);
            double molecule_prob = std::stod(splitted[5]);
            double x1 = std::stod(splitted[6]);
            double y1 = std::stod(splitted[7]);
            double z1 = std::stod(splitted[8]);
            double x2 = std::stod(splitted[9]);
            double y2 = std::stod(splitted[10]);
            double z2 = std::stod(splitted[11]);
            default_mol.push_back({x1, y1, z1});
            default_mol.push_back({x2, y2, z2});
            box.add_molecule_type(elements, rc, molecule_prob, default_mol);
        }
        else if (argc == 13){ // three atoms, two dimensions
            elements.push_back(splitted[3]);
            elements.push_back(splitted[4]);
            double rc = std::stod(splitted[5]);
            double molecule_prob = std::stod(splitted[6]);
            double x1 = std::stod(splitted[7]);
            double y1 = std::stod(splitted[8]);
            double x2 = std::stod(splitted[9]);
            double y2 = std::stod(splitted[10]);
            double x3 = std::stod(splitted[11]);
            double y3 = std::stod(splitted[12]);
            default_mol.push_back({x1, y1});
            default_mol.push_back({x2, y2});
            default_mol.push_back({x3, y3});
            box.add_molecule_type(elements, rc, molecule_prob, default_mol);
        }
        else if (argc == 16){ // three atoms, three dimensions
            elements.push_back(splitted[3]);
            elements.push_back(splitted[4]);
            double rc = std::stod(splitted[5]);
            double molecule_prob = std::stod(splitted[6]);
            double x1 = std::stod(splitted[7]);
            double y1 = std::stod(splitted[8]);
            double z1 = std::stod(splitted[9]);
            double x2 = std::stod(splitted[10]);
            double y2 = std::stod(splitted[11]);
            double z2 = std::stod(splitted[12]);
            double x3 = std::stod(splitted[13]);
            double y3 = std::stod(splitted[14]);
            double z3 = std::stod(splitted[15]);
            default_mol.push_back({x1, y1, z1});
            default_mol.push_back({x2, y2, z2});
            default_mol.push_back({x3, y3, z3});
            box.add_molecule_type(elements, rc, molecule_prob, default_mol);
        }
        else{
            std::cout << "'add moleculetype' does not support " + std::to_string(argc - 3) + " arguments. Aborting." << std::endl;
            exit(0);
        }
    }
    */
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


/* -------------------------------------------------------------
   This function takes care of the 'take snapshot'-command. 
   Thinking about other commands that would make sense to add
   here, but have no idea right now
---------------------------------------------------------------- */

void take(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 1);
    std::string keyword = splitted[1];
    if(keyword == "snapshot"){
        assert(argc > 2);
        std::string filename = splitted[2];
        //box.snapshot(filename);
    }
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


/* ---------------------------------------------------------------
   This function takes care of the run-commands, including 
   run_mc and run_md. 
------------------------------------------------------------------ */

void run(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 2);
    std::string keyword = splitted[1];
    int nsteps = std::stod(splitted[2]);
    if(keyword == "mc"){
        assert(argc > 3);
        int nmoves = std::stod(splitted[3]);
        system.run_mc(nsteps, nmoves);
    }
    /*
    else if(keyword == "md"){
        system.run_md(nsteps);
    }
    */
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


/* ---------------------------------------------------------------
   This function takes care of the thermo output command
------------------------------------------------------------------ */

void thermo(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 3);
    int freq = stoi(splitted[1]);
    std::string filename = splitted[2];
    std::vector<std::string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    //box.set_thermo(freq, filename, outputs);
}


/* -------------------------------------------------------------- 
   This function takes care of the dump output command
----------------------------------------------------------------- */

void dump(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 3);
    int freq = std::stoi(splitted[1]);
    std::string filename = splitted[2];
    std::vector<std::string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    //box.set_dump(freq, filename, outputs);
}


/* -------------------------------------------------------------- 
   This function takes care of the write output command
----------------------------------------------------------------- */

void write(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 2);
    std::string keyword = splitted[1];
    if (keyword == "nsystemsize") {
        std::string filename = splitted[2];
        box.write_nsystemsize(filename);
    }
    else {
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


void special(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert (argc > 1);
    std::string keyword = splitted[0];
    if (keyword == "new") {
        std::string opt = splitted[1];
        if (opt == "box") {
            // set second box as default box
        }
    }
}


/* -------------------------------------------------------------- 
   This function takes care of the remove output command
----------------------------------------------------------------- */

void rm(System &system, Box &box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 1);
    std::cout << "Remove is not yet implemented in the parser" << std::endl;
}
