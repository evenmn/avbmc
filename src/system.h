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

    void set_sampler(class Sampler *);
    void set_sampler(const std::string &, std::function<double(int)>);
    void set_rng(class RandomNumberGenerator *);
    void set_rng(const std::string &);
    void set_forcefield(class ForceField *, int = -1);
    void set_forcefield(const std::string &, const std::vector<std::string> &, int = -1);
    void set_forcefield(const std::string &, const std::string &, int = -1);
    void set_boundary(class Boundary *, int = -1);
    void set_boundary(const std::string &, std::valarray<double> = {}, int = -1);
    void add_constraint(class Constraint *, int = -1);
    void snapshot(const std::string &, int = 0);
    void set_dump(int, const std::string &, const std::vector<std::string> &, int = 0);
    void set_thermo(int, const std::string &, const std::vector<std::string> &, int = 0);
        
    void add_move(class Moves *, double);
    void add_box(class Box *);
    void add_box();

    void check_masses();
    int get_maxiter(int);
    void print_logo();
    void print_info();
    void print_mc_info();

    void run_md(int);
    void run_mc(int, int = 1);


    ~System();

    // variables
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;

    bool logo_printed, rng_allocated_externally, sampler_allocated_externally;
    unsigned int nbox, ndim, nmove, nprocess, step, rank;
    double temp, chempot, poteng, time;
    std::string working_dir;

    std::vector<class Box *> boxes;
    //std::vector<std::string> unique_labels;
    //std::vector<std::string> mass_labels;
    //std::vector<double> masses;
    std::vector<class Moves *> moves;
    std::vector<double> moves_prob;
};
