#include "mersennetwister.h"

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator()
{}

double MersenneTwister::next_gaussian(double mean, double variance)
{
    /* Returns a double drawn from a Gaussian 
     * distribution with mean 'mean' and variance
     * 'variance'
     */
    std::normal_distribution<double> dis(mean, variance);
    return dis(generator);
}

int MersenneTwister::next_int(int upper_limit)
{
    /* Returns an integer between 0 and 'upper_limit'
     */
    std::uniform_int_distribution<int> dis(0, upper_limit - 1);
    return dis(generator);
}

double MersenneTwister::next_double()
{
    /* Returns a double between 0 and 1 drawn 
     * from a uniform distribution
     */
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}


int MersenneTwister::choice(std::vector<double> probabilities)
{
    /* Inspired by numpy.random.choice, but simplified as we always
     * want to draw only one sample. The samples are
     * also always range(len(probabilities), so we don't
     * need to take them as an argument
     */

    // assert that sum of probabilities is 1
    double sum_probs = 0;
    for (double prob : probabilities){
        sum_probs += prob;
    }
    assert(abs(sum_probs - 1.0) < 0.01);

    double r = next_double();
    double p = 0.0;
    for (int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    return 0;
}
