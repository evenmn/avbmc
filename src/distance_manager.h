#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <valarray>


class DistanceManager
{
public:
    DistanceManager(class Box *, double = 1e-2);
    unsigned int add_cutoff(double);
    unsigned int add_cutoff(double, std::string, std::string);
    void clear_neigh(unsigned int);
    void update_neigh(unsigned int, unsigned int, double);
    double normsq(std::valarray<double>);
    void initialize();
    void update_trans(unsigned int);
    void update_remove(unsigned int);
    void update_insert(unsigned int);
    
    std::vector<std::vector<std::vector<int> > > neigh_lists;
    std::vector<std::vector<double> > distance_mat;
    std::vector<std::vector<std::valarray<double> > > distance_cube;

private:
    double cutoff_tol;
    unsigned int ncutoff;
    std::vector<double> cutoffs;
    std::vector<bool> res; // whether or not cutoff is restricted
    std::vector<unsigned int> types1, types2;
    class Box *box = nullptr;
};


