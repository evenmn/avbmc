#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <map>
#include <memory>


class System
{
public:
    System(const std::string & = "."); 
    System(const System &);

    // methods
    void set_sampler(class Sampler*);
    void set_rng(class RandomNumberGenerator*);
        
    void set_temp(double);
    void set_chempot(double);
    void set_mass(std::string, double);

    void add_move(class Moves *, double);
    void add_box(class Box *);

    void check_masses();
    int get_maxiter(int);
    void print_logo();
    void print_info();
    void print_mc_info();

    void run_md(int);
    void run_mc(int, int);

    ~System();

    // variables
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;
    //std::unique_ptr<class RandomNumberGenerator> rng = nullptr;

    bool logo_printed;
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
