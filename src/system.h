#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <memory>
#include <map>


class System
{
public:
    System(std::string  = ""); 

    // methods
    void set_forcefield(class ForceField*);
    //void set_integrator(class Integrator* integrator_in);
    void set_sampler(class Sampler*);
    void set_rng(class RandomNumberGenerator*);
        
    void set_temp(double);
    void set_chempot(double);
    void set_mass(std::string, double);

    //void add_move(class Moves*, double);
    //void add_move(class Moves, double);
    void add_move(std::shared_ptr<class Moves>, double);
    void add_molecule_type(std::string, double);
    void add_molecule_type(std::vector<std::string>, double, double, std::vector<std::valarray<double> >);
    //void add_box(class Box*);
    void add_box(class Box *);
    //void add_box(std::shared_ptr<class Box>);

    void check_masses();
    void init_simulation();
    void init_molecules();
    int get_maxiter(int);
    void print_logo();
    void print_info();
    void print_mc_info();

    void run_md(int);
    void run_mc(int, int);

    ~System();

    // variables
    class ForceField* forcefield = nullptr;
    //class Integrator* integrator = nullptr;
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;
    class MoleculeTypes* molecule_types = nullptr;

    bool initialized;
    int nbox, ndim, ntype, nmove, nprocess, step, rank;
    double temp, chempot, poteng, time;
    std::string working_dir;

    std::map<std::string, int> label2type;
    //std::vector<std::shared_ptr<class Box> > boxes;
    std::vector<class Box *> boxes;
    std::vector<std::string> unique_labels;
    std::vector<std::string> mass_labels;
    std::vector<double> masses;
    std::vector<std::shared_ptr<class Moves> > moves;
    //std::vector<class Moves> moves;
    std::vector<double> moves_prob;
};
