#include "mersennetwister.h"

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator()
{}

//Mersenne Twister RNG
random_device rd;   // will be used to obtain a seed for the random number engine
mt19937 gen(rd());  // standard mersenne_twister_engine seeded with rd()

double MersenneTwister::next_gaussian(double mean, double variance)
{
    /* Returns a double drawn from a Gaussian 
     * distribution with mean 'mean' and variance
     * 'variance'
     */
    normal_distribution<double> dis(mean, variance);
    return dis(gen);
}

int MersenneTwister::next_int(int upper_limit)
{
    /* Returns an integer between 0 and 'upper_limit'
     */
    uniform_int_distribution<int> dis(0, upper_limit - 1);
    return dis(gen);
}

double MersenneTwister::next_double()
{
    /* Returns a double between 0 and 1 drawn 
     * from a uniform distribution
     */
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}

/*
auto MersenneTwister::choice(vector<typename iterator_traits<Iter>::value_type> samples, vector<double> probabilities)
{
    // Inspired by numpy.random.choice
    //
{
    // assert that sum of probabilities is 1
    double sum_probs = 0;
    for (double prob : probabilities){
        sum_probs += prob;
    }
    assert((sum_probs - 1.0) < 0.01);
    assert(samples.size() == probabilities.size());

    // pick sample based on probabilities
    double r = next_double();
    double p = 0.0;
    for (int i=0; i<samples.size(); i++) {
        p += probabilities[i];
        if (r>=p) return samples[i];
    }
}
*/

int MersenneTwister::choice(vector<double> probabilities)
{
    /* Inspired by numpy.random.choice, but simplified as we always
     * want to draw only one sample. The samples are
     * also always range(len(probabilities), so we don't
     * need to take them as an argument
     */

    double r = next_double();
    double p = 0.0;
    for (int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r>=p) return i;
    }
}
