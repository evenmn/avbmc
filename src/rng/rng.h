#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <armadillo>
#include <cassert>

using namespace std;
using namespace arma;

class RandomNumberGenerator
{
public:
    RandomNumberGenerator() : generator(seed()) {}
    virtual int next_int(int upper_limit) = 0;
    virtual double next_double() = 0;
    virtual double next_gaussian(double mean=0., double variance=1.) = 0;
    virtual int choice(vector<double> probabilities) = 0;
    //virtual auto choice(vector<typename iterator_traits<Iter>::value_type> samples, vector<double> probabilities) = 0;
    //virtual mat random_uniform_matrix(int nrows, int ncols) = 0;
    //virtual mat random_normal_matrix(Eigen::Index row, Eigen::Index col, double mean, double variance) = 0;

    virtual ~RandomNumberGenerator() = default;

    random_device seed;
    mt19937 generator;
};
