#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int next_int(int upper_limit);
    double next_double();
    double next_gaussian(double mean, double variance);
    int choice(vector<double> probabilities);
    //auto choice(vector<typename iterator_traits<Iter>::value_type> samples, vector<double> probabilities);
    //Eigen::MatrixXd randomUniformMatrix(Eigen::Index row, Eigen::Index col);
    //Eigen::MatrixXd randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance);
};
