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
#include <map>
#include <functional>


class System
{
public:
    System(const std::string & = ".", bool = true); 
    System(const System &);

    // methods
    void set_temp(double);
    void set_chempot(double);
    void set_mass(std::string, double);
    void set_seed(unsigned int);

    void set_sampler(class Sampler *);
    void set_sampler(const std::string &, std::function<double(int)>, int = 100);
    void set_sampler(const std::string &, std::valarray<double> = {});
    void set_rng(class RandomNumberGenerator *);
    void set_rng(const std::string &);
    void set_forcefield(class ForceField *, int = -1);
    void set_forcefield(const std::string &, const std::vector<std::string> &, int = -1);
    void set_forcefield(const std::string &, const std::string &, int = -1);
    void set_boundary(class Boundary *, int = -1);
    void set_boundary(const std::string &, std::valarray<double> = {}, int = -1);
    void add_constraint(class Constraint *, int = -1);
    void add_constraint(const std::string &, const std::string &, const std::string &, double, int = 0, int = 0);
    void set_dump(int, const std::string &, const std::vector<std::string> &, int = 0);
    void set_thermo(int, const std::string &, const std::vector<std::string> &, int = 0);
    void snapshot(const std::string &, int = 0);
        
    void add_move(class Moves *, double);
    void add_move(const std::string &, double = 1., double = 0.1, double = 0.1, const std::string & = "", int = -1); // for trans and transmh
    void add_move(const std::string &, double, const std::string &, double = 0.95, double = 3.0, bool = false, int = 0, int = 1); // for avbmc, avbmcin, avbmcout, avbmcswap, avbmcswapright
    void add_move(const std::string &, double, std::vector<class Particle>, double = 1.3, double = 0.95, double = 3.0, bool = false, bool = false, int = 0, int = 1); // for avbmcmol, avbmcmolin, avbmcmolout, avbmcmolswap, avbmcmolswapright
    void add_box(class Box *);
    void add_box();
    void add_particle(class Particle, int = 0);
    void add_particle(const std::string &, std::valarray<double>, int = 0);
    void add_particles(std::vector<class Particle>, int = 0);
    void add_particles(const std::string &, std::vector<std::valarray<double> >, int = 0);
    void read_particles(const std::string &, int = 0);
    void write_size_histogram(const std::string &, int = 0);
    std::vector<unsigned int> get_size_histogram(int = 0);

    void rm_particle(unsigned int, int = 0);
    void rm_box(unsigned int);
    void rm_move(unsigned int);
    void rm_constraint(unsigned int, int = 0);
    void clear_particles(int = -1);

    void check_masses();
    unsigned int get_maxiter(unsigned int);
    void print_logo();
    void print_info();
    void print_mc_info();

    void run_md(unsigned int);
    void run_mc(unsigned int, unsigned int = 1);
    void initialize_mc_run();
    void run_mc_cycle(unsigned int = 1);


    ~System();

    // variables
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;

    bool logo_printed, rng_allocated_externally, sampler_allocated_externally;
    int nbox;  // nbox is signed as it is compared to box_id which is signed
    unsigned int ndim, nmove, step;
    double temp, chempot, poteng, time;
    std::string working_dir;

    std::vector<class Box *> boxes;
    //std::vector<std::string> unique_labels;
    //std::vector<std::string> mass_labels;
    //std::vector<double> masses;
    std::vector<class Moves *> moves;
    std::vector<double> moves_prob;
    std::vector<bool> moves_allocated_in_system;
};
